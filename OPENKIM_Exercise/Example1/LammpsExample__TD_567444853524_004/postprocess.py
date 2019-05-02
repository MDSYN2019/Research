"""
Author: Sang Young Noh 
Version: 0.0.1 
Date: 27-04-2019 (Date is in British format so feel free to change if necessary!)
Title: Develop a Python Tool for Generating KIM Property Instances from LAMMPS Output

Description:

All aspects for the KIM system are designed to be extended by user contributions.
However, there are a number of aspects of the system that are currently challenging 
for potential contributors. 

In principle a KIM Test can be as simple as writing a LAMMPS input script to compute a material 
property, such as the face-centered cubic (FCC) equilibrium lattice constant. However, in pracrise 
such a test must also include a script that will post-process the LAMMPS output file to extract 
all necessary dataand put it into a form that the rest of the KIM system can readily unserstand. 

Although there are lots of material scientists who could create and contribute interesting simulations,
few are capable of creating the necessary post-processing software to complete a contribution to KIM.
Thus, it is crucial for the success for KIM that we create a suite of software tools that make it easy 
for the contributors to extract heir from LAMMPS output files and import it into he apporpriate KIM format.
 

My Personal Notes:

The LAMMPS version used in this version is "lammps-12Dec18"

=> Key point - How could we generalize this python processing code to work for a wider variety of LAMMPS input files?

1. Its probably best to divide into blocks - at the moment, the lammps input file is divided into multiple blocks:
   
   - Units and atom_style block, maybe add boundary. Could name this 'conditions' 
   - Box properties and dimensions - 
   - Output style - dump (pdb, dcd, ...) 
   - Variable definition
   - Print definitions

These have to be divided into blocks and made into modules so that they can be puzzled together, with error condiitons making sure 
incompatible fixes are rooted out.

2.  The program needs to be modular - not complicated enough that it is just boilerplate in the end and also easy for even novices 
     in Python to update for their own purposes and for future versions of lammps
"""

import argparse
import string
import os
import re  # regex module 
import unittest
import subprocess
import signal
import psutil, time
import edn_format
import pprint
from lammps import lammps

LammpsPath = "/home/oohnohnoh1/Desktop/LAMMPS/lammps-12Dec18/src"
KIMModelsPath = "/usr/local/lib/kim-api/models"
TemplatePath = os.path.abspath('../')
CurrentPath = os.path.abspath('.')

KIMMODELSLIST = [
	'LennardJones612_UniversalShifted__MO_959249795837_003',
	'LennardJones_Ar',
	'ex_model_Ar_P_LJ',
	'SW_StillingerWeber_1985_Si__MO_405512056662_005',
	'ex_model_Ar_P_MLJ_Fortran',
	'ex_model_Ar_P_Morse',
	'ex_model_Ar_P_Morse_07C',
	'ex_model_Ar_P_Morse_07C_w_Extensions',
	'ex_model_Ar_P_Morse_MultiCutoff',
	'ex_model_Ar_SLJ_MultiCutoff'
]

PROPERTY = [
	'atomic-mass',
	'bulk-modulus-isothermal-cubic-crystal-npt',
	'bulk-modulus-isothermal-hexagonal-crystal-npt',
	'cohesive-energy-lattice-invariant-shear-path-cubic-crystal',
	'cohesive-energy-lattice-invariant-shear-unrelaxed-path-cubic-crystal',
	'cohesive-energy-relation-cubic-crystal',
	'cohesive-energy-shear-stress-path-cubic-crystal',
	'cohesive-free-energy-cubic-crystal',
	'cohesive-free-energy-hexagonal-crystal',
	'cohesive-potential-energy-2d-hexagonal-crystal',
	'cohesive-potential-energy-cubic-crystal',
	'cohesive-potential-energy-hexagonal-crystal',
	'configuration-cluster-fixed',
	'configuration-cluster-relaxed',
	'configuration-nonorthogonal-periodic-3d-cell-fixed-particles-fixed',
	'configuration-nonorthogonal-periodic-3d-cell-fixed-particles-relaxed',
	'configuration-nonorthogonal-periodic-3d-cell-relaxed-particles-fixed',
	'configuration-nonorthogonal-periodic-3d-cell-relaxed-particles-relaxed',
	'configuration-periodic-2d-cell-fixed-particles-fixed',
	'elastic-constants-first-strain-gradient-isothermal-cubic-crystal-npt',
	'elastic-constants-first-strain-gradient-isothermal-monoatomic-hexagonal-crystal-npt',
	'elastic-constants-isothermal-cubic-crystal-npt',
	'enthalpy-of-mixing-curve-substitutional-binary-cubic-crystal-npt',
	'enthalpy-of-mixing-curve-substitutional-binary-cubic-crystal-nvt',
	'extrinsic-stacking-fault-relaxed-energy-fcc-crystal-npt',
	'gamma-surface-relaxed-fcc-crystal-npt',
	'grain-boundary-symmetric-tilt-energy-ideal-cubic-crystal',
	'grain-boundary-symmetric-tilt-energy-relaxed-cubic-crystal',
	'grain-boundary-symmetric-tilt-energy-relaxed-relation-cubic-crystal',
	'intrinsic-stacking-fault-relaxed-energy-fcc-crystal-npt',
	'linear-thermal-expansion-coefficient-cubic-crystal-npt',
	'melting-temperature-constant-pressure-cubic-crystal',
	'monovacancy-formation-energy-monoatomic-cubic-diamond',
	'monovacancy-neutral-formation-free-energy-crystal-npt',
	'monovacancy-neutral-migration-energy-crystal-npt',
	'monovacancy-neutral-relaxation-volume-crystal-npt',
	'monovacancy-neutral-relaxed-formation-potential-energy-crystal-npt',
	'monovacancy-neutral-unrelaxed-formation-potential-energy-crystal-npt',
	'monovacancy-neutral-unrelaxed-formation-potential-energy-crystal-npt',
	'phonon-dispersion-dos-cubic-crystal-npt',
	'phonon-dispersion-relation-cubic-crystal-npt',
	'shear-stress-path-cubic-crystal',
	'stacking-fault-relaxed-energy-curve-fcc-crystal-npt',
	'structure-2d-hexagonal-crystal-npt',
	'structure-cubic-crystal-npt',
	'structure-hexagonal-crystal-npt',
	'structure-monoclinic-crystal-npt',
	'structure-orthorhombic-crystal-npt',
	'structure-rhombohedral-crystal-npt',
	'structure-tetragonal-crystal-npt',
	'structure-triclinic-crystal-npt',
	'surface-energy-broken-bond-fit-cubic-bravais-crystal-npt',
	'surface-energy-cubic-crystal-npt',
	'surface-energy-ideal-cubic-crystal',
	'unstable-stacking-fault-relaxed-energy-fcc-crystal-npt',
	'unstable-twinning-fault-relaxed-energy-fcc-crystal-npt',
	'verification-check'
]

PROPERTYDENVAL = ['instance-id',
 'short-name',
 'species',
 'a',
 'basis-atom-coordinates',
 'space-group',
 'wyckoff-multiplicity-and-letter',
 'wyckoff-species',
 'wyckoff-coordinates',
 'cohesive-potential-energy']

#os.mkdir(str(TemplatePath + "/" + "output")) # Make output directory 
# Parser to read in the force fields, the Lattice constant, log file and the lammps binary of the Computer

#parser = argparse.ArgumentParser(description="Stdin for OpenKIM")
#parser.add_argument('--Forcefield', action = 'store' , type = str, help = 'The name of the KIM forcefield')
#parser.add_argument('--Lattice_Constant', action = 'store', type = float, help = 'The value of the lattice constant')
#parser.add_argument('--Log', action = 'store', type = str, help = 'Name of the log file')
#parser.add_argument('--LAMMPS_binary', action='store', type = str, help = 'Path to lammps binary')
#print (parser.parse_args(['--Forcefield', '--Lattice_Constant', '--Log', '--LAMMPS_binary']))

# Making sure of errors

# Make sure that the input Forcefield exists

#try:
#	args.Forcefield in KIMMODELSLIST
#except IndexError:
#	print ("Input forcefield is not valid")

# Testing that LAMMPS works through a subprocess

#try:
#	proc = subprocess.Popen("/home/oohnohnoh1/Desktop/LAMMPS/lammps-12Dec18/src/lmp_ubuntu", shell=True) 
#	p_id = psutil.Process(proc.pid)
#	p_id.kill()
#	print ("LAMMPS works!")
#except subprocess.CalledProcessError as e:
#	print(e.output)

# Testing whether log.lammps output is in the output folder 

#try:
#	FileOpen = open(CurrentPath + "/" + "log.lammps")
#	print ("log.lammps exists")
#except FileNotFoundError:
#	print ("log.lammps does not exist")

"""
Example for the property definition of the cohesive energy relation of a cubic crystal

{
  "property-id" "tag:staff@noreply.openkim.org,2014-04-15:property/cohesive-energy-relation-cubic-crystal"

  "property-title" "Cohesive energy versus lattice constant relation for a cubic crystal"

  "property-description" "Cohesive energy versus lattice constant relation for a cubic crystal at zero absolute temperature.  Lattice constants are taken to correspond to the conventional cubic unit cell.  Moreover, note that here the cohesive energy is defined as the *negative* of the potential energy per atom."

  "short-name" {
    "type"         "string"
    "has-unit"     false
    "extent"       [":"]
    "required"     false
    "description"  "Short name defining the cubic crystal type."
  }
  "species" {
    "type"         "string"
    "has-unit"     false
    "extent"       [":"]
    "required"     true
    "description"  "The element symbols of the basis atoms.  The order in which the species are specified must correspond to the order of the atoms listed in 'basis-atom-coordinates'."
  }
  "a" {
    "type"         "float"
    "has-unit"     true
    "extent"       [":"]
    "required"     true
    "description"  "A vector of conventional unit cell lattice constants of the cubic crystal."
  }
  "basis-atom-coordinates" {
    "type"         "float"
    "has-unit"     false
    "extent"       [":",3]
    "required"     true
    "description"  "Fractional coordinates of the basis atoms in the conventional unit cell.  If the unit cell vectors are denoted by <a>, <b>, and <c>, and the fractional coordinates of atom 'i' are [afrac_i, bfrac_i, cfrac_i], the value of 'basis-atom-coordinates' will be of the form [[afrac_1 bfrac_1 cfrac_1] [afrac_2 bfrac_2 cfrac_2] ... ].  All components of each basis atom should be between zero and one, inclusive of zero."
  }
  "space-group" {
    "type"         "string"
    "has-unit"     false
    "extent"       []
    required       false
    "description"  "Hermann-Mauguin designation for the space group associated with the symmetry of the crystal (e.g. Immm, Fm-3m, P6_3/mmc)."
  }
  "wyckoff-multiplicity-and-letter" {
    "type"         "string"
    "has-unit"     false
    "extent"       [":"]
    required       false
    "description"  "Multiplicity and standard letter of Wyckoff sites (e.g. 4a, 2b).  The order of elements in this array must correspond to the order of the entries listed in 'wyckoff-species' and 'wyckoff-coordinates'."
  }
  "wyckoff-species" {
    "type"         "string"
    "has-unit"     false
    "extent"       [":"]
    required       false
    "description"  "The element symbol of the atomic species of the unique Wyckoff sites used in the fully symmetry-reduced description of the crystal.  The order of elements in this array must correspond to the order of the entries listed in 'wyckoff-multiplicity-and-letter' and 'wyckoff-coordinates'."
  }
  "wyckoff-coordinates" {
    "type"         "float"
    "has-unit"     false
    "extent"       [":",3]
    required       false
    "description"  "Coordinates of the Wyckoff sites, given as fractions of the lattice vectors.  The order of elements in this array must correspond to the order of the entries listed in 'wyckoff-species' and 'wyckoff-multiplicity-and-letter'."
  }
  "cohesive-potential-energy" {
    "type"         "float"
    "has-unit"     true
    "extent"       [":"]
    "required"     true
    "description"  "Cohesive energy (negative of the potential energy per atom) associated with the corresponding lattice constant."
  }
"""
		
 GenericName = ["short-name", "species", "a", "basis-atom-coordinates", "space-group",  "wyckoff-multiplicity-and-letter", "wyckoff-species", "wyckoff-coordinates", "cohesive-potential-energy"]

 GenericProperties = ["type", "has-unit", "extent", required", "description"]


def EdnSourceValue(key):
	"""
	The KIM infrastructure embraces a subset of EDN as a standard data format. EDN stands for extensible data notation, and is 
	pronounced like the word "eden"	
	"""
	# First, check if the property is indeed a valid one as defined in OpenKIM (https://openkim.org/properties)
	try:
		key in PROPERTY
	except KeyError:
		print ("Input property not found")

	propertyArray = ["short-name", "species", "a", "basis-atom-coordinates", "space-group", "wyckoff-multiplicity-and-letter", "wyckoff-species", "wyckoff-coordinates", "cohesive-potential-energy"] 

	# Output dictionary part
	
	outputDict = {}
	outputDict["property-id"] = "tag:staff@noreply.openkim.org,2014-04-15:property/{}".format(key)	
	for property in propertyArray:
		outputDict[property] = {}
		outputDict[property]["source-value"] = None
		outputDict[property]["source-unit"] = None

    ## An example output following the case for the result.edn.tpl as shown in the OpenKIM example

	# short-name
	outputDict["short-name"]["source-value"] = ["fcc"]

	# "species" 
	outputDict["species"]["source-value"] = ["Ar", "Ar", "Ar", "Ar"]

	# "a"
	outputDict["a"]["source-value"] = 3.60000000000017
	outputDict["a"]["source-unit"] = "angstrom"

	# basis-atom-coordinates
	outputDict["basis-atom-coordinates"]["source-value"] = [[0, 0, 0],
															 [0, 0.5, 0.5],
															 [0.5, 0, 0.5],
															 [0.5, 0.5, 0]
	                                                          ]
	# space-group
	outputDict["space-group"]["source-unit"] = "Fm-3m"

	# wyckoff-species
	outputDict["wyckoff-species"]["source-unit"] = ["Ar"]

	# wyckoff-potential
	outputDict["wyckoff-coordinates"]["source-value"] = [[0,0,0]]

    # cohesive-potential-energy
	outputDict["cohesive-potential-energy"]["source-value"] = -6.46911040131582
	outputDict["cohesive-potential-energy"]["source-unit"] =  "eV"

	# Check for entries that have not been used, and remove it if it contains 'None'
	for property in propertyArray:
		outputDict[property] = {k: v for k, v in outputDict[property].items() if v is not None}

	# Convert output to a edn format
	ednOutput = edn_format.dumps(outputDict)
	return ednOutput

def printEdn(output):
	with open("example.out", "w+") as fout:
		fout.write(output)			

class KIMPostprocess:
	"""

	Class to read in lammps.log and template files and call to produce the edn files	

	"""
	def __init__(self, logfile, input_template, writefile, path):
		"""
		Constructor 
		"""
		try:
			self.LogfileInput = open(str(path + "/" + logfile), "rb")
			self.InputTemplateInput = open(str(path + "/" + input_template), "rb") # Need to rename this 
			self.LogfileRead = self.logfile_input.readlines()
			self.InputTemplateInput = self.input_template_input.readlines() # Need to rename this 
		except IOError:
			print ("Make sure you have the right paths for the lammps and template files")
	def PropertySearch(self):
		"""
		"""
		with open("../lammps.in.template","r+") as fin:
			filedata = fin.read()
			filedata = filedata.replace("sed_initial_lattice_constant_string", args.Lattice_Constant)
			filedata = filedata.replace("sed_model_string", args.Forcefield)
		with open("lammps.in", "w+") as fout:
			fout.write(filedata)			
	def EdnWriter(self):		
		"""
		Writer for the Edn file 
		"""
		FinalPressureLine = [line for line in self.logfile_read.split(' ') if "Final Pressure" in line] 
		CohesiveEnergyLine = [line for line in self.logfile_read.split(' ') if "Cohesive Energy" in line]
		LatticeConstantLine = [line for line in self.logfile_read.split(' ') if "lattice constant" in line]

		self.FinalPressureVal = FinalPressureLine[0].decode('utf-8').rstrip().split('=')
		self.CohesiveEnergyVal = CohesiveEnergyLine[0].decode('utf-8').rstrip().split('=')
		self.LatticeConstantVal = LatticeConstantLine[0].decode('utf-8').rstrip().split('=')

		with open("../results.edn.tpl","r+") as fin:
			filedata = fin.read()
			filedata = filedata.replace("_LATCONST_", self.LatticeConstantVal[0])
			filedata = filedata.replace("_ECOHESIVE_", self.CohesiveEnergyVal[0])
		with open("results.edn", "w+") as fout:
			fout.write(filedata)
	
