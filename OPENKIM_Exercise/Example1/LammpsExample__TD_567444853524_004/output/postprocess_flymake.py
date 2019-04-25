"""

Sang Young Noh - OpenKIM exercise for lammps - For Professor Tadmor 

"""
import argparse
import string
import os
import re  # regex module 
KIM_MODELS_DIR = "/usr/local/lib/kim-api/models"


# Instead of the sed echo commands, we can use argparse

current_path = os.path.abspath('.')
os.mkdir(str(current_path + "/" + output)) # Make output directory 

parser = argparse.ArgumentParser(description='Stdin for OpenKIM')
parser.add_argument('Forcefield', metavar = 'FF' , type = str, help = 'The name of the KIM forcefield')
parser.add_argument('Lattice_Constant', metavar = 'C', type = float, help = 'The value of the lattice constant')
parser.add_argument('Log', metavar = 'L', type = str, help = 'Name of the log file')
parser.add_argument('lammps_binary', metavar = 'B', type = str, help = 'Path to lammps binary')
args = parser.parse_args()
print(args.accumulate(args.integers))

# Regular expression patterns

I_pattern = re.compile("sed_initial_lattice_constant_string")
II_pattern = re.compile("<(\d{4,5})>")

s = open("/home/oohnohnoh1/Desktop/GIT/MD_Design_and_Research/OPENKIM_Exercise/Example1/LammpsExample__TD_567444853524_004/lammps.in.template","r+")
for i, line in enumerate(s.readlines()):
	if re.search(I_pattern, line):
		print ("Found on line {}: {}".format(i, line))

class KIM_Postprocess:
	def __init__(self, logfile, input_template,path):
		self.logfile_input = open(str(path + "/" + logfile))
		self.input_template_input = open(str(path + "/" + input_template)) # Need to rename this 
		self.logfile_read = self.logfile_input.readlines()
		self.input_template_input = self.input_template_input.readlines() # Need to rename this 
	def propertySearch(self):
		model_string_pattern = re.compile("sed_model_string")
		lattice_contant_pattern = re.compile("sed_initial_lattice_constant_string")
		# Extract values
		finalpressure_line = [line for line in self.logfile_read.split(' ') if "Final Pressure" in line] 
		ecohesive_line = [line for line in self.logfile_read.split(' ') if "Cohesive Energy" in line]
		latticeconstant_line = [line for line in self.logfile_read.split(' ') if "lattice constant" in line]
	def logfileReader(self):
		pass
	def output(self):
		pass

	
# How to generalize the LAMMPS input file?

"""
1. Its probably best to divide into blocks - at the moment, the lammps input file is divided into multiple blocks:
   
   - Units and atom_style block, maybe add boundary. Could name this 'conditions' 

   - Box properties and dimensions - 
 
   - Output style - dump (pdb, dcd, ...) 
   
   - Variable definition

   - Print definitions

These have to be divided into blocks and made into modules so that they can be puzzled together, with error condiitons making sure 
incompatible fixes are rooted out.

"""
