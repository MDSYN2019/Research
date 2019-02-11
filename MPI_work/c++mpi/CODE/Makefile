#
# GENERAL C++ MAKEFILE
#

.KEEP_STATE:
.SUFFIXES: .cpp .H .h

# Change EXEC to the name of the executable which you wish to
# create. Make OFILES be the name of all the .C files in your program,
#

ifdef MPISRC
LIB    =  libSCmathlibP.a
else
LIB    =  libSCmathlib.a
endif

OFILES	= SCmathlib    \
	SCchapter2     \
	SCchapter3     \
	SCchapter4     \
	SCchapter5     \
	SCchapter6     \
	SCchapter7     \
	SCchapter8     \
	SCchapter9     \
	SCchapter10


PROGRAMS = \
 
#Compiler specifications
CCX = g++
ifdef MPISRC
DEBUG  =  -O3 -DMPISRC 
else
DEBUG  =  -g
endif
IFLAGS	= -I.
LFLAGS	= -L. 
CFLAGS  = $(DEBUG) 
LIBS	=  -lm 
RUN_FLAGS = 	
AR = ar rv
LDFLAGS = $(LFLAGS) -lSC $(LIBS)

# These are different make options.  You can invoke them by typing 
# the label name.  For example to time your make, type 'make time' -- 
# 'make' being the name of this command and 'time' being the label used 
# below to define the result.  
#

default: all

all: lib $(PROGRAMS)

lib: $(LIB)

.cpp.o:
	$(CCX) $(CFLAGS) -c $<

LMOD  = $(foreach routine, $(OFILES), $(routine).o)

$(LIB): $(LMOD)
	$(AR) $(LIB) $(LMOD)

# Clean up after making. 
clean:
	$(RM) *.o $(LIB) $(PROGRAMS)






