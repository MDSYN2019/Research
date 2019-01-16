FCOMPFLAGS    =	-O
CCOMPFLAGS    =	-O

FFLAGS        =	$(FCOMPFLAGS)
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	f77
FC            =	f77
CC            =	cc

MAKEFILE      =	Makefile
PROGRAM       =	barrier1

INCS	      =	system.inc

OBJS	      =	force.o             \
		integrate_and.o     \
		integrate_mc.o      \
		integrate_nve.o     \
		integrate_res.o     \
		barrier1.o          \
		mdloop.o            \
		ran_uniform.o       \
		ran_vel.o           \
		ran_gauss.o         \
		readdat.o           \
		sample_diff.o       \
		sample_prof.o       \
		sstmm.o

all:		$(PROGRAM)

$(PROGRAM)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROGRAM)::	$(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJS) -o $(PROGRAM)

clean:;		@rm -f $(OBJS) core
