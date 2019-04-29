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
 

Notes:

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

# Optional and currently experiementing - running lammps on python

#from lammps import lammps
#lmp = lammps()
#lmp.file("in.lj")

LAMMPSPATH = "/home/oohnohnoh1/Desktop/LAMMPS/lammps-12Dec18/src"
KIMMODELSDIR = "/usr/local/lib/kim-api/models"
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
	
TemplatePath = os.path.abspath('../')
CurrentPath = os.path.abspath('.')
os.mkdir(str(template_path + "/" + "output")) # Make output directory 

# Parser to read in the force fields, the Lattice constant, log file and the lammps binary of the Computer

parser = argparse.ArgumentParser(description="Stdin for OpenKIM")
parser.add_argument('--Forcefield', action = 'store' , type = str, help = 'The name of the KIM forcefield')
parser.add_argument('--Lattice_Constant', action = 'store', type = float, help = 'The value of the lattice constant')
parser.add_argument('--Log', action = 'store', type = str, help = 'Name of the log file')
parser.add_argument('--LAMMPS_binary', action='store', type = str, help = 'Path to lammps binary')
args = parser.parse_args()
print (args)

# Making sure of errors

try:
	args.Forcefield in KIM_MODELS_LIST
except IndexError:
	print ("Input forcefield is not valid")

# Could try to use the python version of LAMMPS here for exceptions 

try:
	proc = subprocess.Popen("/home/oohnohnoh1/Desktop/LAMMPS/lammps-12Dec18/src/lmp_ubuntu", shell=True) 
	p_id = psutil.Process(proc.pid)
	p_id.wait(timeout = 3)
except subprocess.CalledProcessError as e:
	print(e.output)

class KIM_Postprocess:
	"""

	Class to read in lammps.log and template files and call to produce the edn files	

	"""
	def __init__(self, logfile, input_template, writefile, path):
		"""
		Constructor 
		"""
		self.LogfileInput = open(str(path + "/" + logfile), "rb")
		self.InputTemplateInput = open(str(path + "/" + input_template), "rb") # Need to rename this 
		self.LogfileRead = self.logfile_input.readlines()
		self.InputTemplateInput = self.input_template_input.readlines() # Need to rename this 
	def PropertySearch(self):
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

			
# How to generalize the LAMMPS input file?


""" 

Test block - using pytest 

"""
def test_input():
	pass

def test_output():
	pass
