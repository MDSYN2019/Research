FCOMPFLAGS    =	-O
CCOMPFLAGS    =	-O

FFLAGS        =	$(FCOMPFLAGS)
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	f77
FC            =	f77
CC            =	cc

MAKEFILE      =	Makefile
PROGRAM       =	md_nve

INCS	      =	system.inc

OBJS	      =	force.o             \
		init.o              \
		integrate.o         \
		md.o                \
		mdloop.o            \
		ran_uniform.o       \
		rangauss.o          \
		readdat.o           \
		sample_diff.o       \
		sample_gyra.o       \
		sstmm.o             \
		writepdb.o

all:		$(PROGRAM)

$(PROGRAM)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROGRAM)::	$(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJS) -o $(PROGRAM)

clean:;		@rm -f $(OBJS) core
