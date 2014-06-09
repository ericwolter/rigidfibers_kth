TARGET = ADVECT_FIBERS
HOME =/Users/eric/Developer/kth/rigidfibers/fortran

#  -------------- PUT IN DIRECTORIES HERE THAT YOU WANT ------
BINDIR = $(HOME)

# Location of source files (could be your working directory)

FORLIBS.f = $(HOME)
FORLIBS.f90 = $(HOME)
LIBDIR = $(HOME)
LIBDIR2=/opt/NAG/fnsol03d9
INCLS=/opt/NAG/fnsol03d9/nag_mod_dir
MOD_DIR = $(HOME)
RM       = rm -f
MV       = mv -f

F90FLAGS += -g -Og
#	    -C\
           # -xtarget=ultra2i -xcache=16/32/1:2048/64/1 \
	   # -xlibmopt \
	   # -O3

F77FLAGS += -O3


LAPACKDIR = /usr/lib
LOADLIBES = -L$(LAPACKDIR) -llapack -lblas

COMPILE.f90= gfortran $(F90FLAGS) -I. -I$(BINDIR) 
LINK.f90   = gfortran $(F90FLAGS)
COMPILE.f  = gfortran $(F77FLAGS)
LINK       = gfortran $(F90FLAGS)

#  Files we don't have to build:
.PHONY:  makefile depends.mak  $(TARGET).objs


#  The dependencies for $(TARGET) are in a special file:
include $(TARGET).objs

#
#  Pattern rules
#
vpath %.o $(BINDIR)
vpath %.f $(FORLIBS.f)
vpath %.f90 . $(FORLIBS.f90) $(MOD_DIR)

%.o:: %.f
	$(COMPILE.f) -c $< -o $@
	@echo -----------

%.o:: %.f90
	$(COMPILE.f90) -c $< -o $@
	@echo -----------

%.mod:: %.f90
	$(COMPILE.f90) -c $< -o $@
	@echo -----------



$(TARGET) : $($(TARGET).objs)
	$(LINK) $(filter %.o, $^) $(LOADLIBES) -o $@
#	-$(MV) *.o $(BINDIR)

cleanobj:
	-$(RM) core $(BINDIR)/*.o *.mod *.o

#
#  Make a tags file for regexp replacements and searching
#
TAGS: $(wildcard *.f90) $(wildcard *.fi) $(wildcard *.objs)
	etags *.f90 *.fi *.objs

# Dependencies files
#include depends.mak
