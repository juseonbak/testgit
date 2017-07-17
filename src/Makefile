##########################################################
#                                                         #
#   This is the Makefile for SOMIPROF, the SAO's          #
#                                                         #
#     BOas REtrieval for Atmospheric Spectra              #
#                                                         #
#                                                         #
#                                                         #
#   Author:  Thomas P. Kurosu, Kelly Chance               #
#                                                         #
#            Harvard-Smithsonian Center for Astrophysics  #
#            Atomic and Molecular Physics Division        #
#            60 Garden Street (MS 50)                     #
#            Cambridge, MA 02138, USA                     #
#            Tel/Fax: +1 -- 617 - 495 7213                #
#            EMail tkurosu@cfa.harvard.edu                #
#                                                         #
#                                                         #
#   Last modified: June 2003                              #
#                                                         #
###########################################################

# ================================================
# Define some environment variables that
# are unlikely to change between host-types
# ================================================
AR      = ar
ARFLAGS = ruv
MV      = mv
RM      = rm -f
SHELL   = bash

# HDF/HDF-EOS/Toolkit/Pspline libraries
PSPLINELIBS  = -lpspline -lezcdf
NETCDFLIBS   = -lnetcdf
HE4LIBS      = -lhdfeos -lmfhdf -lGctp -ldf -ljpeg -lsz -lz
HE5LIBS      = -lhdf5_fortran -lhe5_hdfeos -lhdf5 -lhdf5hl_fortran -lhdf5_hl 
PGSLIBS      = -lPGSTK
OMILIBS      = -lOMIutil_v2.1.0

#===============================================================


# =======================================================
# Target Directories for Modules, Library, and Executable
# =======================================================
PGETMPDIR  = $(PGEHOME)/tmp
PGEBINDIR  = $(PGEHOME)/bin
PGEOUTDIR  = $(PGEIODIR)
PGEOBJDIR  = $(PGEHOME)/obj

OMSAOLIB   = $(SAOPGE)
MAINDIR    = omi-pge
MAIN       = $(OMSAOLIB)_main
TARGET     = $(SAOPGE).exe
EXECUTABLE = $(TARGET) 

#PGSHOME = $(OMIUTIL)
#PGSBIN  = $(OMIUTIL)/bin/$(FLAVOR)

include ../make/make.flags

# =============================================================
#  Library for all compiled object files. Rather than keeping
#  the *.o in the current directory, they are archived in
#  LIBFILE. This way, we can maintain object files for various
#  compilation flags (opt,dbx,...) or even architectures
# =============================================================
LOCALLIB = $(PGETMPDIR)/lib$(OMSAOLIB).a
VPATH = ./$(PGEBINDIR):$(PGETMPDIR):$(PGEHOME)/src:$(PGEHOME)/src/gome:\
	$(PGEHOME)/src/omi-pge:$(PGEHOME)/src/o3prof:$(PGEHOME)/src/omps-pge

# ================================================
# All sources are "sourced out" to an include file
# ================================================
include ../make/make.sources
include ../make/make.gome_sources
#include ../make/make.omps_sources
include ../make/make.o3prof_sources


# =================================================
#  Everything we want to make in this distribution
# =================================================
all : envcheck dirs messages mainlink $(EXECUTABLE)
	@for target in $^; do                          \
           $(MAKE) $$target;     \
        done

# ============================================
#  Here's how the111133.executable is made from the
#  archive (LOCALLIB) and the various Modules
# ============================================
$(EXECUTABLE) : $(LOCALLIB) $(MAIN).f90 $(F90SOURCES) $(F90GOMESOURCES) $(F90OMPSOURCES) $(OZSOURCES) \
                $(LIDORTPATH)/$(LIDORTLIBFILE) $(PELIBPATH)/$(PELIBFILE)
	$(FC) $(MAIN).f90 -o $(TARGET) $(LNFLAGS) $(FFLAGS) $(IFLAGS) $(LFLAGS) $(LIBS)$ # $(GEMS_LIB)


# ======================================
#  The SAO_OMOCLO object library archive
# ======================================
$(LOCALLIB): $(OBJECTS) $(GOMEOBJECTS)  $(OZOBJECTS) $(OMPSOBJECTS)


# =============================================================
# A list of directories that are empty at CVS checkin, and are
# thus not version tagged (files only). They must be created at
# the start of the PGE make process.
# =============================================================
emptydirs = $(PGETMPDIR) $(PGSMSG) $(PGEBINDIR) $(PGEOUTDIR)

# =====================
#  Various Make targets
# =====================
build: all

install: all
	$(MAKE) cleanlog cleanoutput
	cp -p $(EXECUTABLE) $(PGEOUTDIR)/
	mv 	*.mod $(PGETMPDIR)


test: all
	$(MAKE) install
	cd $(PGEOUTDIR); $(PGEOUTDIR)/$(EXECUTABLE) 2> $(SAOPGE).stderr > $(SAOPGE).stdout
	@grep -i "exit code" $(PGEOUTDIR)/LogStatus

vtest: all
	$(MAKE) install
	cd $(PGEOUTDIR); $(PGEOUTDIR)/$(EXECUTABLE)
	@grep -i "exit code" $(PGEOUTDIR)/LogStatus

run:
	cd $(PGEOUTDIR); $(PGEOUTDIR)/$(EXECUTABLE)

include ../make/make.targets
include ../make/make.patterns
#include ../make/make.omps_patterns
include ../make/make.o3prof_patterns
include ../make/make.gome_patterns

# =================
#  End of Makefile
# =================
