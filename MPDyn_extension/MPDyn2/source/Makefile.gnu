########################################
##
## written by W. Shinoda
##
########################################

## 5) for Intel compiler ver 11.0
FC     = gfortran
CC     = gcc
#$FFLAGS = -w -static -O3 -tpp7 -ipo -tune pn4 -pad -cpp -DPCC -DPUBFFT -DCUBIC -DSCINT
FFLAGS = -ffree-form -Wall -O -x f95-cpp-input -DPCC -DPUBFFT -DCUBIC -DSCINT -DCOMPRESS -DSURF
LFLAGS = 
PFLAGS = -ffree-form -Wall -O -x f95-cpp-input -DPCC -DPUBFFT -DCUBIC -DSCINT
#SFLAGA = -L/opt/intel/mkl61/lib/32/ -lmkl_lapack -lmkl_ia32 -lguide -lpthread
#PFLAGA = -L/opt/intel/mkl61/lib/32/ -lmkl_lapack -lmkl_ia32 -lguide -lpthread
SFLAGA = 
PFLAGA = 
PFC    = mpif90

include make.mpdyn

#< end of Makefile>
