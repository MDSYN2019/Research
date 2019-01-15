FCOMPFLAGS    =	-O
CCOMPFLAGS    =	-O

FFLAGS        =	$(FCOMPFLAGS)
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	f77
FC            =	f77
CC            =	cc

MAKEFILE      =	Makefile
PROGRAM       =	mc_nvt

INCS	      =	system.inc conf.inc parameter.inc potential.inc

OBJS	      =	adjust.o     \
		ener.o       \
		eneri.o      \
		lattice.o    \
		mc_nvt.o     \
		mcmove.o     \
		readdat.o    \
		sample.o     \
		sstmm.o      \
		store.o      \
		toterg.o     \
		writepdb.o   \
		ran_uniform.o

all:		$(PROGRAM)

$(PROGRAM)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROGRAM)::	$(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJS) -o $(PROGRAM)

clean:;		@rm -f $(OBJS) core
