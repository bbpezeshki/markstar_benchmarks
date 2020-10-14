import tools
import osprey
osprey.start()

# load the confspace
confspaces = tools.load_confspaces(
    '../designs/disc_binding_ss/design_00000.toml')
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
parallelism = osprey.Parallelism(cpuCores=4)
ecalc = osprey.EnergyCalculator(complexConfSpace, ffparams, parallelism=parallelism)

# MARK* needs a rigid energy calculator too
ecalcRigid = osprey.SharedEnergyCalculator(ecalc, isMinimizing=False)

# configure K*
kstar = osprey.KStar(
    proteinConfSpace,
    ligandConfSpace,
    complexConfSpace,
    epsilon=0.683, 
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
scoredSequences = kstar.run(ecalc.tasks)

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
