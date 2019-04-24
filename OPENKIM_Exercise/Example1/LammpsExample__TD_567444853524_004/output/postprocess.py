"""

Sang Young Noh - OpenKIM exercise for lammps

"""
import argparse

parser = argparse.ArgumentParser(description='Stdin for OpenKIM')
parser.add_argument('Forcefield', metavar = 'FF' , type = str, help = 'The name of the KIM forcefield')
parser.add_argument('Lattice Constant', metavar = 'C', type = float, help = 'The value of the lattice constant')
args = parser.parse_args()
print(args.accumulate(args.integers))

class KIM_Postprocess:
	def __init__(self):
		pass
	
# Do I need inheritance or not..?
