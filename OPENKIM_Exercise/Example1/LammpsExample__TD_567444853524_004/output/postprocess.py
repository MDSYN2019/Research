"""
Author: Sang Young Noh -
Version: 0.0.1 
Date: --::--::----
Title: Develop a Python Tool for Generating KIM Property Instances from LAMMPS Output

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


KIM_MODELS_DIR = "/usr/local/lib/kim-api/models"

# Temporary placeholder 
KIM_MODELS_LIST = [
'LennardJones612_UniversalShifted__MO_959249795837_003',
'LennardJones_Ar',
'SW_StillingerWeber_1985_Si__MO_405512056662_005',
'ex_model_Ar_P_LJ',
'ex_model_Ar_P_MLJ_Fortran',
'ex_model_Ar_P_Morse',
'ex_model_Ar_P_Morse_07C',
'ex_model_Ar_P_Morse_07C_w_Extensions',
'ex_model_Ar_P_Morse_MultiCutoff',
'ex_model_Ar_SLJ_MultiCutoff'
]
	
template_path = os.path.abspath('../')
current_path = os.path.abspath('.')


os.mkdir(str(current_path + "/" + "output")) # Make output directory 

# Instead of the sed echo commands, we can use argparse

parser = argparse.ArgumentParser(description="Stdin for OpenKIM")
parser.add_argument('--Forcefield', action = 'store' , type = str, help = 'The name of the KIM forcefield')
parser.add_argument('--Lattice_Constant', action = 'store', type = float, help = 'The value of the lattice constant')
parser.add_argument('--Log', action = 'store', type = str, help = 'Name of the log file')
parser.add_argument('--LAMMPS_binary', action='store', type = str, help = 'Path to lammps binary')
args = parser.parse_args()
print (args)

try:
	args.Forcefield in KIM_MODELS_LIST
except IndexError:
	print ("Input forcefield is not valid")



s = open("/home/oohnohnoh1/Desktop/GIT/MD_Design_and_Research/OPENKIM_Exercise/Example1/LammpsExample__TD_567444853524_004/lammps.in.template","r+")
for i, line in enumerate(s.readlines()):
	if re.search(I_pattern, line):
		print ("Found on line {}: {}".format(i, line))

class KIM_Postprocess:

	def __init__(self, logfile, input_template, writefile, path):

		self.logfile_input = open(str(path + "/" + logfile), "rb")

		self.input_template_input = open(str(path + "/" + input_template), "rb") # Need to rename this 

		self.logfile_read = self.logfile_input.readlines()

		self.input_template_input = self.input_template_input.readlines() # Need to rename this 

	def propertySearch(self):

		model_string_pattern = re.compile("sed_model_string")

		lattice_contant_pattern = re.compile("sed_initial_lattice_constant_string")

		# Extract values

		finalpressure_line = [line for line in self.logfile_read.split(' ') if "Final Pressure" in line] 

		ecohesive_line = [line for line in self.logfile_read.split(' ') if "Cohesive Energy" in line]

		latticeconstant_line = [line for line in self.logfile_read.split(' ') if "lattice constant" in line]

	def edn_writer(self):

		self.writefile = writefile

		f = open(str(self.writefile), "wb")

		# do something

		f.close()

	def output(self):

		pass
# How to generalize the LAMMPS input file?


""" 

Test block - using pytest 

"""
def test_input():
	pass

def test_output():
	pass
