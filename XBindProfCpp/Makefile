#FLAGS = -lm -L/usr/pub/lib -lefence -o
#CC = cc -O2 -fullwarn -TENV:large_GOT
#CC = cc -g -Wall
CC = g++ -O2    # at least for SunOS
#CC = cc -g

#CC = cc -O2 -fullwarn

#CC = cc -O2
FLAGS = -lm  -static

SOURCE = main.cpp interfaceMSA/IfAlignment.cpp interfaceMSA/IfMutant.cpp interfaceMSA/ResName.cpp math/SubstitutionMatrix.cpp utils/StringTools.cpp interfaceMSA/IfEnergyCalculator.cpp  interfaceMSA/IfSingleMutant.cpp  math/AAProbabilityArray.cpp  utils/FileTools.cpp

OBJECT = ${SOURCE:.cpp=.o} 

BINDIR = .

.c.o: 
	$(CC) -c $< -o $@


xbindprof : $(OBJECT)
	$(CC) $(OBJECT) $(FLAGS)  -o $(BINDIR)/xbindprof${ARCH}



clean:
	rm -f $(OBJECT) $(BINDIR)/xbindprof${ARCH}

show:
	echo $(SOURCE)

