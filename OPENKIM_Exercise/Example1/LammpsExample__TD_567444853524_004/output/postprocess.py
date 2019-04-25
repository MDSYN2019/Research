"""

Sang Young Noh - OpenKIM exercise for lammps - For Professor Tadmor 

"""
import argparse
import string
import os
import re  # regex module 


# Instead of the sed echo commands, we can use argparse

current_path = os.path.abspath('.')
os.mkdir(str(current_path + "/" + output)) # Make output directory 

#parser = argparse.ArgumentParser(description='Stdin for OpenKIM')
#parser.add_argument('Forcefield', metavar = 'FF' , type = str, help = 'The name of the KIM forcefield')
#parser.add_argument('Lattice_Constant', metavar = 'C', type = float, help = 'The value of the lattice constant')
#parser.add_argument('Log', metavar = 'L', type = str, help = 'Name of the log file')
#parser.add_argument('lammps_binary', metavar = 'B', type = str, help = 'Path to lammps binary')

#args = parser.parse_args()
#print(args.accumulate(args.integers))

# Regular expression patterns

I_pattern = re.compile("sed_initial_lattice_constant_string")
II_pattern = re.compile("<(\d{4,5})>")

s = open("/home/oohnohnoh1/Desktop/GIT/MD_Design_and_Research/OPENKIM_Exercise/Example1/LammpsExample__TD_567444853524_004/lammps.in.template","r+")
for i, line in enumerate(s.readlines()):
	if re.search(I_pattern, line):
		print ("Found on line {}: {}".format(i, line))
# Replacing sed model string



class KIM_Postprocess:
	def __init__(self, logfile, input_template,path):
		self.logfile = open(str(path + "/" + logfile))
		self.input_template = open(str(path + "/" + input_template))
		
	def logfileReader(self):
		pass
	def output(self):
		pass
	def propertySearch(self):
		finalpressure = re.search()
		ecohesive = re.search()
		latticeconstant = re.search()


# How to generalize the LAMMPS input file?

"""
1. Its probably best to divide into blocks - 


"""
