#include <stdlib.h>
#include "interfaceMSA/ResName.h"
#include "math/SubstitutionMatrix.h"
#include "math/AAProbabilityArray.h"
#include "interfaceMSA/IfSingleMutant.h"
#include "interfaceMSA/IfMutant.h"
#include "interfaceMSA/IfAlignment.h"
#include "interfaceMSA/IfEnergyCalculator.h"
#include "utils/FileTools.h"
#include <iostream>
#include <fstream>

int main(int argc, char** argv)
{
    if (argc!=4) {
        printf("xbindprof align.out mutList 0.5\n");
        return argc;
    }
	string alignFile = argv[1];
	string mutFile = argv[2];
	double cutoff = atof(argv[3]);

	IfAlignment align = IfAlignment(alignFile,cutoff);
	fstream f(argv[2],ios::in);
	if (!f)
	{
		cerr << "can't open file: " << mutFile << endl;
	}
	string line;
	
	while (!f.eof())
	{
		getline(f, line);
                if(line.length() < 4)
                   continue;
		IfMutant mut = IfMutant(line);
		double ddG = 0;
		for (int i = 0; i < mut.mutNum; i++)
		{
			IfSingleMutant smut = mut.mutList.at(i);
			ddG += IfEnergyCalculator::getDDG(align, smut);
		}
		std::printf("%-7.3f\n", ddG);
	}
	f.close();
        return 0;
}
