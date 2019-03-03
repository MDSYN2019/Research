# MD_Jarsynski

This is the C programming script I used to calculate the free energy change of aa nanoparticle translocation, through a model bilayer.
The code uses the work distribution from the trajectory, over 50 runs from a position of 40 \AA (the centre/normal of the bilayer is located 
at 0 \AA) and trajects towards -40 \AA and that is it's final position. All you need to get the script working is to compile the cmake 
file, usually with the following procedure: 

mkdir build
cmake ../. 
make

After, you should get a executable file named 'JE'. I've included a example pdb trajectory of a hydrophobic nanoparticle to demonstrate 
what kind system is being simulated here. 

