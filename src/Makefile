TARGET = ADVECT_FIBERS
HOME =/afs/nada.kth.se/home/c2m2/katarina/Desktop

#  -------------- PUT IN DIRECTORIES HERE THAT YOU WANT ------
BINDIR = $(HOME)/RIGIDFIBERS_NON_PERIODIC_CODE

# Location of source files (could be your working directory)

FORLIBS.f = $(HOME)/RIGIDFIBERS_NON_PERIODIC_CODE
FORLIBS.f90 = $(HOME)/RIGIDFIBERS_NON_PERIODIC_CODE
LIBDIR = $(HOME)/RIGIDFIBERS_NON_PERIODIC_CODE
LIBDIR2=/opt/NAG/fnsol03d9
INCLS=/opt/NAG/fnsol03d9/nag_mod_dir
MOD_DIR = $(HOME)/RIGIDFIBERS_NON_PERIODIC_CODE
RM       = rm -f
MV       = mv -f

F90FLAGS += -O3 -fns\
#	    -C\
           # -xtarget=ultra2i -xcache=16/32/1:2048/64/1 \
	   # -xlibmopt \
	   # -O3

F77FLAGS += -O3


LAPACKDIR = /usr/lib
LOADLIBES = -L$(LAPACKDIR) -llapack -lblas

COMPILE.f90= f90 $(F90FLAGS) -I. -I$(BINDIR) 
LINK.f90   = f90 $(F90FLAGS)
COMPILE.f  = f90 $(F77FLAGS)
LINK       = f90 $(F90FLAGS)

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
include depends.mak
