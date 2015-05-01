

CC       = h5pcc
#CCFLAGS  = -I/opt/local/include -O2 #-ffast-math
CCFLAGS  = -Wall -std=c99 -I/ope/local/include -Ofast -ftree-vectorizer-verbose=1 -fopenmp
#CCFLAGS  = -g -Wall -std=c99 -O0 -save-temps
#CCFLAGS+= -I/usr/local/Cellar/gsl/1.16/include
#CCFLAGS  = -I/opt/local/include -pg

EXTRALIBS = -lm

CC_COMPILE  = $(CC) $(CCFLAGS) -c 
CC_LOAD     = $(CC) $(CCFLAGS) 

.c.o:
	$(CC_COMPILE) $*.c

EXE = harm
all: $(EXE) 

SRCS = \
mpi.c bounds.c coord.c diag.c fixup.c \
init.c interp.c main.c metric.c \
phys.c ranc.c step_ch.c io.c xdmf_output.c \
advance_particles.c init_particles.c pdump.c \
synchrotron.c utoprim_mm.c current_calc.c wind.c

OBJS = \
mpi.o bounds.o coord.o diag.o fixup.o \
init.o interp.o main.o metric.o \
phys.o ranc.o step_ch.o io.o xdmf_output.o \
advance_particles.o init_particles.o pdump.o \
synchrotron.o utoprim_mm.o current_calc.o wind.o

INCS = decs.h  defs.h  

$(OBJS) : $(INCS) makefile

$(EXE): $(OBJS) $(INCS) makefile
	$(CC_LOAD) $(OBJS) $(EXTRALIBS) -o $(EXE)

clean:
	/bin/rm -f *.o 

newrun:
	/bin/rm -rf dumps images/*.ppm ener.out
	/bin/rm -f *.o
	make harm

