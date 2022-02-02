from pathlib import Path
import argparse
import time

import _submodules.tools as tools
import osprey
osprey.start()


arg_parser = argparse.ArgumentParser();
arg_parser.add_argument("--omit-wt-amino-acids", action='store_true', help="signal to not have OSPREY automatically add wt amino acids");
arg_parser.add_argument("--epsilon", type=float, default=0.683, help="bbk* epsilon to use");
arg_parser.add_argument("--num-solutions", type=int, default=1, help="number of solutions for bbk* to compute");
arg_parser.add_argument("--num-cores", type=int, default=1, help="number of CPU cores to use")


args = arg_parser.parse_args();

print("")

if args.omit_wt_amino_acids:
	ADD_WT_AMINO_ACIDS = False
else:
	ADD_WT_AMINO_ACIDS = True
print(">>>>> ADD_WT_AMINO_ACIDS:", str(ADD_WT_AMINO_ACIDS))

EPSILON = args.epsilon; #default: 0.683;
print(">>>>> EPSILON:", str(EPSILON))

NUM_SOLUTIONS = args.num_solutions
print(">>>>> NUM_SOLUTIONS:", str(NUM_SOLUTIONS))

NUM_CORES = args.num_cores
print(">>>>> NUM_CORES:", str(NUM_CORES))

print("\n=============================================")
print("=============================================\n")





##############################################################
def purgeLeftOverOspreyFiles():
	p_purge = Path().absolute();

	SUFFIXES_TO_KEEP = {".py", ".txt", ".stdout"}

	for g_purge in p_purge.glob('*'):
		if g_purge.is_file():
			if g_purge.suffix not in SUFFIXES_TO_KEEP:
				g_purge.unlink();
###############################################################

purgeLeftOverOspreyFiles();

for design_file in Path('designs_to_run_on/').glob('*'):
	if design_file.suffix != ".toml":
		continue;

	design_id = design_file.stem.split("_",1)[1]
	pdb = None;
	with design_file.open('r') as fin:
		for line in fin:
			if line.startswith("molecule"):
				pdb = Path().absolute().parent/"resources"/"pdb" / line.split("=")[1].strip().replace('"',"")
	assert(pdb)
	print()
	print("DESIGN: " + design_id)
	print("PDB: " + pdb.name)
	print()

	# load the confspace
	# confspaces = tools.load_confspaces(
	#     '../designs/disc_binding_ms/design_00000.toml')
	confspaces = tools.load_confspaces(design_file, force_wt=ADD_WT_AMINO_ACIDS)

	ffparams = confspaces['ffparams']

	# read a PDB file for molecular info
	# make sure all strands share the same template library
	templateLib = osprey.TemplateLibrary(ffparams.forcefld)

	# make the conf space for the protein
	proteinConfSpace = confspaces['protein']

	# make the conf space for the ligand
	ligandConfSpace = confspaces['ligand']

	# make the conf space for the protein+ligand complex
	complexConfSpace = confspaces['complex']

	# how should we compute energies of molecules?
	# (give the complex conf space to the ecalc since it knows about all the templates and degrees of freedom)
	parallelism = osprey.Parallelism(cpuCores=NUM_CORES)
	minimizingEcalc = osprey.EnergyCalculator(complexConfSpace, ffparams, parallelism=parallelism, isMinimizing=True)

	# BBK* needs a rigid energy calculator too, for multi-sequence bounds on K*
	rigidEcalc = osprey.SharedEnergyCalculator(minimizingEcalc, isMinimizing=False)

	# how should we define energies of conformations?
	def confEcalcFactory(confSpace, ecalc):
		eref = osprey.ReferenceEnergies(confSpace, ecalc)
		return osprey.ConfEnergyCalculator(confSpace, ecalc, referenceEnergies=eref)

	# configure BBK*
	bbkstar = osprey.BBKStar(
		proteinConfSpace,
		ligandConfSpace,
		complexConfSpace,
		numBestSequences=NUM_SOLUTIONS,
		epsilon=EPSILON, 
		writeSequencesToConsole=True,
		writeSequencesToFile='bbkstar.results.tsv'
	)

	# configure BBK* inputs for each conf space
	for info in bbkstar.confSpaceInfos():

		# how should we define energies of conformations?
		eref = osprey.ReferenceEnergies(info.confSpace, minimizingEcalc)
		info.confEcalcMinimized = osprey.ConfEnergyCalculator(info.confSpace, minimizingEcalc, referenceEnergies=eref)

		# compute the energy matrix
		ematMinimized = osprey.EnergyMatrix(info.confEcalcMinimized, cacheFile='emat.%s.dat' % info.id)

		# how should confs be ordered and searched?
		# (since we're in a loop, need capture variables above by using defaulted arguments)
		def makeAStarMinimized(rcs, emat=ematMinimized):
			return osprey.AStarTraditional(emat, rcs, showProgress=False)
		info.confSearchFactoryMinimized = osprey.BBKStar.ConfSearchFactory(makeAStarMinimized)

		# BBK* needs rigid energies too
		confEcalcRigid = osprey.ConfEnergyCalculatorCopy(info.confEcalcMinimized, rigidEcalc)
		ematRigid = osprey.EnergyMatrix(confEcalcRigid, cacheFile='emat.%s.rigid.dat' % info.id)
		def makeAStarRigid(rcs, emat=ematRigid):
			return osprey.AStarTraditional(emat, rcs, showProgress=False)
		info.confSearchFactoryRigid = osprey.BBKStar.ConfSearchFactory(makeAStarRigid)

		# how should we score each sequence?
		# (since we're in a loop, need capture variables above by using defaulted arguments)
		# def makePfunc(rcs, confEcalcMinimized=info.confEcalcMinimized, ematMinimized=ematMinimized, confEcalcRigid=confEcalcRigid, ematRigid=ematRigid):
		# 	return osprey.MARKStarPfunc(
		# 		confEcalcMinimized.confSpace,
		# 		ematMinimized,
		# 		confEcalcMinimized,
		# 		ematRigid,
		# 		confEcalcRigid,
		# 		rcs
		# 	)
		# info.pfuncFactory = osprey.KStar.PfuncFactory(makePfunc)
		def makePfunc(rcs, confEcalc=info.confEcalcMinimized, emat=ematMinimized):
			return osprey.PartitionFunction(
				confEcalc,
				osprey.AStarTraditional(emat, rcs, showProgress=False),
				osprey.AStarTraditional(emat, rcs, showProgress=False),
				rcs
			)
		info.pfuncFactory = osprey.KStar.PfuncFactory(makePfunc)

	# run BBK*
	# scoredSequences = bbkstar.run(minimizingEcalc.tasks)
	print()
	print("\t----- Calling: bbkstar.run() -----")
	print()
	start = time.time()
	scoredSequences = bbkstar.run(minimizingEcalc.tasks)
	end = time.time()

	# make a sequence analyzer to look at the results
	analyzer = osprey.SequenceAnalyzer(bbkstar)

	# use results
	for scoredSequence in scoredSequences:
		print("result:")
		print("\tsequence: %s" % scoredSequence.sequence)
		print("\tK* score: %s" % scoredSequence.score)

		# write the sequence ensemble, with up to 10 of the lowest-energy conformations
		numConfs = 10
		analysis = analyzer.analyze(scoredSequence.sequence, numConfs)
		print(analysis)
		analysis.writePdb(
			'seq.%s.pdb' % scoredSequence.sequence,
			'Top %d conformations for sequence %s' % (numConfs, scoredSequence.sequence)
		)

	# use results
	for i,scoredSequence in enumerate(scoredSequences):
		print()
		print("result " + str(i+1) + ":")
		print("\tsequence: %s" % scoredSequence.sequence)
		print("\tscore: %s" % scoredSequence.score)
	print()
	print("ELAPSED TIME: " + str(end - start))

	print()
	print("DESIGN: " + design_id)
	print("PDB: " + pdb.name)
	print("\n=============================================")
	print("=============================================\n")

	purgeLeftOverOspreyFiles();

purgeLeftOverOspreyFiles();