TARGET = fibers

# The compiler
FC = gfortran

# https://gcc.gnu.org/ml/fortran/2009-05/msg00409.html
# Debug flags
#FCFLAGS = -g -fbounds-check -ftree-vectorizer-verbose=2 -ffree-form -ffree-line-length-none
# Release flags
# similar options as ifort defaults
FCFLAGS = -O3 -ffast-math -funroll-loops --param max-unroll-times=4 -mtune=native -ffree-form -ffree-line-length-none

# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
#LDFLAGS =

# List of executables to be built within the package
PROGRAMS = fibers

# "make" builds all
all: $(PROGRAMS)

fibers.o: constants.incl

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.F; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.fpp
	$(FC) $(FCFLAGS) -c $<

%.o: %.FPP
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)
