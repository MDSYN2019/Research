FCOMPFLAGS    =	-O
CCOMPFLAGS    =	-O

FFLAGS        =	$(FCOMPFLAGS)
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	f77
FC            =	f77
CC            = cc

MAKEFILE      =	Makefile
PROG1         =	mc-npt

INCS	      =	system.inc

OBJ1	      =	mc-npt.o writepdb.o ran_uniform.o sstmm.o

all:		$(PROG1)

$(PROG1)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROG1)::	$(OBJ1) $(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJ1) -o $(PROG1)

clean:;		@rm -f $(OBJ1) core
