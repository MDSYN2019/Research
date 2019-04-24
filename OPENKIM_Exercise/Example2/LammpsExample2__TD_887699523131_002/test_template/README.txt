Author: Daniel S. Karls (karl0100 |AT| umn DOT edu), University of Minnesota
Date: 10/10/2014

This example LAMMPS Test uses the LAMMPSExample2__TD_887699523131_002 Test Driver to compute
the cohesive energy of a {{lattice}} lattice of {{fullsymbol}} atoms at a variety of lattice spacings
ranging from {{aminfrac}}*a_0 to {{amaxfrac}}*a_0.
The value of "a_0" (as it is referred to in the Test Driver) passed to the Test Driver is taken to be
the equilibrium lattice constant for the Model being used, as obtained by querying the KIM database for
Test Results from the LatticeConstantCubicEnergy_{{lattice}}_{{symbol}} Test when paired against the relevant Model.

Note: Any Test which makes use of a Test Driver can have as its executable the standard python
script used as the Test executable here.
