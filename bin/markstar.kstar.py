from pathlib import Path
import argparse
import time

import _submodules.tools as tools
import osprey
osprey.start()


arg_parser = argparse.ArgumentParser();
arg_parser.add_argument("--omit-wt-amino-acids", action='store_true', help="signal to not have OSPREY automatically add wt amino acids");
arg_parser.add_argument("--epsilon", type=float, default=0.683, help="k* epsilon to use");
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

    # choose a forcefield
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
    ecalc = osprey.EnergyCalculator(complexConfSpace, ffparams, parallelism=parallelism)

    # MARK* needs a rigid energy calculator too
    ecalcRigid = osprey.SharedEnergyCalculator(ecalc, isMinimizing=False)

    # configure K*
    kstar = osprey.KStar(
        proteinConfSpace,
        ligandConfSpace,
        complexConfSpace,
        epsilon=EPSILON, 
        writeSequencesToConsole=True,
        writeSequencesToFile='kstar.results.tsv'
    )

    # configure K* inputs for each conf space
    for info in kstar.confSpaceInfos():

        # how should we define energies of conformations?
        eref = osprey.ReferenceEnergies(info.confSpace, ecalc)
        info.confEcalc = osprey.ConfEnergyCalculator(info.confSpace, ecalc, referenceEnergies=eref)

        # compute the minimized energy matrix
        ematMinimized = osprey.EnergyMatrix(info.confEcalc, cacheFile='emat.%s.dat' % info.id)

        # if you want to use MARK*, pass in a rigid energy calculator and energy matrix as well
        confEcalcRigid = osprey.ConfEnergyCalculatorCopy(info.confEcalc, ecalcRigid)
        ematRigid = osprey.EnergyMatrix(confEcalcRigid, cacheFile='emat.%s.rigid.dat' % info.id)

        # how should we score each sequence?
        # (since we're in a loop, need capture variables above by using defaulted arguments)
        def makePfunc(rcs, confEcalcMinimized=info.confEcalc, ematMinimized=ematMinimized, confEcalcRigid=confEcalcRigid, ematRigid=ematRigid):
            return osprey.MARKStarPfunc(
                confEcalcMinimized.confSpace,
                ematMinimized,
                confEcalcMinimized,
                ematRigid,
                confEcalcRigid,
                rcs
            )
        info.pfuncFactory = osprey.KStar.PfuncFactory(makePfunc)

    # run K*
    # scoredSequences = kstar.run(ecalc.tasks)
    print()
    print("\t----- Calling: kstar.run() -----")
    print()
    start = time.time()
    scoredSequences = kstar.run(ecalc.tasks)
    end = time.time()

    # make a sequence analyzer to look at the results
    analyzer = osprey.SequenceAnalyzer(kstar)

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

