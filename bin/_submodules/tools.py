#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""tools

Tools related to loading designs from TOML formatted configuration files. Most
methods require a running OSPREY instance. This instance can be started with:
    `osprey.start()`
"""

##################################################
# Version
#   0.0
# Author
#   Graham Holt
#
# Contact Info:
#    Bruce Donald
#    Duke University
#    Department of Computer Science
#    Levine Science Research Center (LSRC)
#    Durham
#    NC 27708-0129
#    USA
#    e-mail: www.cs.duke.edu/brd/
#
# <signature of Bruce Donald>, Mar 1, 2018
# Bruce Donald, Professor of Computer Science
#
##################################################

##################################################
# Imports
##################################################
import osprey
import toml
# from os.path import join, dirname, abspath
from pathlib import Path
import multiprocessing
##################################################
# Globals
##################################################
# app_path = dirname(abspath(__file__))
# PDB = join(app_path, "resources/pdb")
# DESIGN = join(app_path, "designs")
# TEST = join(app_path, "tests")
parent_directory = Path().absolute().parent;
PDB = parent_directory / 'resources/pdb/';

##################################################
# Classes / Functions
##################################################
def load_confspaces(f, xtal_rotamers=True, continuous=False, force_wt=True):
    """Loads OSPREY ConfSpace objects from a TOML file.

        OPTIONS:
            xtal_rotamers:  if true, add the rotamers from the input pdb
            continuous:     if true, set continuous rotamers
            force_wt:       if true, always add the wild-type amino acid

    """
    config = toml.load(f)
    return loadd_confspaces(config,
                           xtal_rotamers=xtal_rotamers,
                           continuous=continuous,
                           force_wt=force_wt)

def loadd_confspaces(config, xtal_rotamers=True, continuous=False, force_wt=True):
    """Loads OSPREY ConfSpace objects from a config dictionary.

        OPTIONS:
            xtal_rotamers:  if true, add the rotamers from the input pdb
            continuous:     if true, set continuous rotamers
            force_wt:       if true, always add the wild-type amino acid

    """
    ffparams = osprey.ForcefieldParams()
    mol = osprey.readPdb(PDB / config["molecule"])
    template_library = osprey.TemplateLibrary(ffparams.forcefld)
    # Make sure we don't have 3 strands (some cfs files do)
    if len(config["strand_definitions"]) > 2:
        exit()


    # Define the protein strand
    protein = osprey.Strand(mol,
                            templateLib=template_library,
                            residues=config["strand_definitions"]["strand0"]
                           )
    for resi, res_allowed in config["strand_mutations"]["strand0"].items():

        # Add the osprey.WILD_TYPE object to the allowed residues if desired
        if force_wt:
            res_allowed.append(osprey.WILD_TYPE)

        # Set the flexibility
        protein.flexibility[resi]\
                .setLibraryRotamers(*res_allowed)
        if xtal_rotamers and force_wt:
            protein.flexibility[resi].addWildTypeRotamers()
        if continuous:
            protein.flexibility[resi].setContinuous()

    # Define the ligand strand
    ligand = osprey.Strand(mol,
                           templateLib=template_library,
                           residues=config["strand_definitions"]["strand1"]
                          )
    for resi, res_allowed in config["strand_mutations"]["strand1"].items():

        # Add the osprey.WILD_TYPE object to the allowed residues if desired
        if force_wt:
            res_allowed.append(osprey.WILD_TYPE)

        # Set the flexibility
        ligand.flexibility[resi]\
                .setLibraryRotamers(*res_allowed)
        if xtal_rotamers and force_wt:
            ligand.flexibility[resi].addWildTypeRotamers()
        if continuous:
            ligand.flexibility[resi].setContinuous()

    # Build spaces
    return {'protein' : osprey.ConfSpace(protein),
            'ligand' : osprey.ConfSpace(ligand),
            'complex': osprey.ConfSpace([protein, ligand]),
            'ffparams' : ffparams
           }

def make_emat(fi, fo):
    """Save complex energy matrix based on TOML file"""
    confspaces = load_confspaces(fi)
    parallelism = osprey.Parallelism(cpuCores=multiprocessing.cpu_count())
    makeo_emat(confspaces["complex"],
               confspaces["ffparams"],
               parallelism,
               fo)

def maked_emat(config, fo):
    """Save complex energy matrix based on dictionary"""
    confspaces = loadd_confspaces(config)
    parallelism = osprey.Parallelism(cpuCores=multiprocessing.cpu_count())
    makeo_emat(confspaces["complex"],
               confspaces["ffparams"],
               parallelism,
               fo)

def makeo_emat(confspace, ffparams, parallelism, fo, continuous=False):
    """Save energy matrix based on OSPREY objects"""
    ecalc = osprey.EnergyCalculator(confspace,
                                    ffparams,
                                    parallelism,
                                    isMinimizing=continuous)
    eref = osprey.ReferenceEnergies(confspace, ecalc)
    conf_ecalc = osprey.ConfEnergyCalculator(confspace,
                                             ecalc,
                                             referenceEnergies=eref)
    emat = osprey.EnergyMatrix(conf_ecalc,
                               cacheFile=fo)
