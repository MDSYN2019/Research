import argparse
import string
import os
import re  # regex module 



parser = argparse.ArgumentParser(description="Stdin for OpenKIM")
parser.add_argument('--Forcefield', action = "store" , type = str, help = 'The name of the KIM forcefield')
parser.add_argument('--Lattice_Constant', action = "store",  type = float, help = 'The value of the lattice constant')
parser.add_argument('--Log', type = str, action = "store",  help = 'Name of the log file')
parser.add_argument('--lammps_binary', action = "store",  type = str, help = 'Path to lammps binary')
args = parser.parse_args()
print (args)

# Need to make sure the Forcei

print ("We are processing the argument inputs with the following: {} ".format(args.Forcefield))

