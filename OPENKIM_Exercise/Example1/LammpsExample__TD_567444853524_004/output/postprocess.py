"""

Sang Young Noh - OpenKIM exercise for lammps

"""
import argparse
import string
 
parser = argparse.ArgumentParser(description='Stdin for OpenKIM')
parser.add_argument('Forcefield', metavar = 'FF' , type = str, help = 'The name of the KIM forcefield')
parser.add_argument('Lattice Constant', metavar = 'C', type = float, help = 'The value of the lattice constant')
parser.add_argument('Log', metavar = 'L', type = str, help = 'Name of the log file')

args = parser.parse_args()
print(args.accumulate(args.integers))

s = open("/home/oohnohnoh1/Desktop/GIT/MD_Design_and_Research/OPENKIM_Exercise/Example1/LammpsExample__TD_567444853524_004/lammps.in.template","r+")
for line in s.readlines():
	print (line)
	
class KIM_Postprocess:
	def __init__(self):
		pass
	
# Do I need inheritance or not..?
