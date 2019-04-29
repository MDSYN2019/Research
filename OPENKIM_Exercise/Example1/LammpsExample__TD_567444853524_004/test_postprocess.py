from lammps import lammps # Importing LAMMPS modules 
import argparse
import string
import os
import re  # regex module 
import unittest
import subprocess
import signal
import psutil, time
import pytest
from postprocess import KIMPostprocess

@pytest.mark.file
def test_something():
	file = open()
	pass

@pytest.mark.somethingelse
def test_something_another():
	pass

