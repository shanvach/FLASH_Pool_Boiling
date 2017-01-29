# FLASH makefile definitions for Cray XT5 System, Kraken at NICS
#
# The XT5 makes use of modules.  
# Load the following modules before compiling (names as of 5th Nov 08).
#
# HDF5: module load hdf5
# PNETCDF: module load pnetcdf 
# TAU: module load tau


#----------------------------------------------------------------------------
# Set the HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

export cur-dir := $(shell pwd)

# Set the location of top directory
export setup_dir = $(cur-dir)

MPI_PATH   =
PAPI_PATH  = 
PAPI_FLAGS = 
NCMPI_PATH = 


#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------


FCOMP   =  ftn
CCOMP   =  cc
CPPCOMP =  CC
LINK    =  ftn


#-----------------------------------------------------------------------------
# Compilation flags
#
#  Three sets of compilation/linking flags are defined: one for optimized code
#  code ("-opt"), one for debugging ("-debug"), and one for testing ("-test").
#  Passing these flags to the setup script will cause the value associated with
#  the corresponding keys (i.e. those ending in "_OPT", "_DEBUG", or "_TEST") to
#  be incorporated into the final Makefile. For example, passing "-opt" to the
#  setup script will cause the flags following "FFLAGS_OPT" to be assigned to
#  "FFLAGS" in the final Makefile. If none of these flags are passed, the default
#  behavior will match that of the "-opt" flag.
#  In general, "-opt" is meant to optimize compilation and linking. "-debug"
#  should enable runtime bounds checking, debugger symbols, and other compiler-
#  specific debugging options. "-test" is useful for testing different
#  combinations of compiler flags particular to your individual system.
#----------------------------------------------------------------------------

# Test Compiler PGI:
#FFLAGS_TEST = -c

# PGI Compiler:
#FFLAGS_OPT   = -c -r8 -fpic -fastsse -Minfo -Mneginfo -tp amd64
#FFLAGS_DEBUG = -c -r8 -g -fpic -Ktrap=divz -Mchkfpstk -Mchkptr -Mchkstk -Mdclchk -Mbounds
#FFLAGS_TEST  = -c -r8 -fpic -Mprof=lines

# GNU Compiler:
#FFLAGS_OPT =  -c -O2 -fdefault-real-8 -fdefault-double-8 \
#-ffree-line-length-none -Wuninitialized
#FFLAGS_DEBUG = -ggdb -c -fdefault-real-8 -fdefault-double-8 \
#-ffree-line-length-none -pedantic -Wall -Wextra -Waliasing \
#-Wsurprising -Wconversion -Wunderflow \
#-ffpe-trap=invalid,zero,overflow -fbounds-check \
#-fbacktrace -fdump-core -finit-real=nan \
#-finit-integer=-999999 -fimplicit-none
#FFLAGS_TEST =  -c -fdefault-real-8 -fdefault-double-8 \
#-ffree-line-length-none

#-ffree-line-length-none -Wuninitialized -optTauSelectFile="/nics/b/home/mvanella/ADAPTATIVE/FLASH4_INS/sites/kraken.nics.tennessee.edu/select.tau"

# CRAY Compiler:
FLAGS_OPT   = -c -O3 -hfp3 -rm -xomp -s real64 -s integer32  \
###-optTauSelectFile="/nics/d/home/kdelane2/gradWork/flash4_mar2012_cleanHYPRE/INS/sites/kraken.nics.tennessee.edu/select.tau"
FLAGS_DEBUG = -c -G 0 -s real64 -s integer32
FLAGS_TEST  = -c  -s real64 -s integer32


# INTEL compiler:
##FFLAGS_OPT   = -c -r8 -i4 -O3 -real_size 64 
#\
#-optTauSelectFile="/nics/b/home/mvanella/ADAPTATIVE/FLASH4_INS/sites/kraken.nics.tennessee.edu/select.tau"
##FFLAGS_DEBUG = -c -g -r8 -i4 -real_size 64
#-O0 -check bounds -check format \
#-check output_conversion -warn all -warn error -real_size 64 -check uninit \
#-traceback -fp-stack-check  -fpe0 -check pointers
##FFLAGS_TEST  = -c -g -r8 -i4


CFLAGS_OPT   = -c -O3
CFLAGS_DEBUG = -c -g
CFLAGS_TEST  = -c

#No path required because we are using compiler wrapper scripts.
FFLAGS_MPI   = 
CFLAGS_MPI   = 

# if we are using HDF5, we need to specify the path to the include files
FFLAGS_HDF5  = -I${CRAY_HDF5_DIR}/hdf5-gnu/include  ${HDF5_FLIB} -DH5_USE_16_API
CFLAGS_HDF5  = -I${CRAY_HDF5_DIR}/hdf5-gnu/include  ${HDF5_CLIB} -DH5_USE_16_API
CFLAGS_NCMPI = ${PNETCDF_LIB}

FFLAGS_PAPI  = 

#CFLAGS_HYPRE = ${HYPRE_LIB}
#FFLAGS_HYPRE = ${HYPRE_LIB}
CFLAGS_HYPRE = -I/sw/xt-cle3.1/hypre/2.7.0b/cnl3.1_gnu4.6.1/include -L/sw/xt-cle3.1/hypre/2.7.0b/cnl3.1_gnu4.6.1/lib -lHYPRE
FFLAGS_HYPRE = -I/sw/xt-cle3.1/hypre/2.7.0b/cnl3.1_gnu4.6.1/include -L/sw/xt-cle3.1/hypre/2.7.0b/cnl3.1_gnu4.6.1/lib -lHYPRE


#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -o  #if using IPA need it on link line too -Mipa=fast,inline
LFLAGS_DEBUG = -g -o
LFLAGS_TEST  = -Mprof=lines -o

#----------------------------------------------------------------------------
# Library specific linking
#
#  If a FLASH module has a 'LIBRARY xxx' line in its Config file, we need to
#  create a macro in this Makefile.h for LIB_xxx, which will be added to the
#  link line when FLASH is built.  This allows us to switch between different
#  (incompatible) libraries.  We also create a _OPT, _DEBUG, and _TEST
#  library macro to add any performance-minded libraries (like fast math),
#  depending on how FLASH was setup.
#----------------------------------------------------------------------------

LIB_HDF5    = -L${CRAY_HDF5_DIR}/hdf5-gnu/lib ${HDF5_CLIB}
#LIB_HDF5    = -L/opt/cray/hdf5/1.8.5.0/hdf5-gnu/lib ${HDF5_CLIB}

              
LIB_MPI     = 
LIB_PAPI    = 
LIB_PNG     = 

LIB_OPT     = 
LIB_DEBUG   =
LIB_TEST    =

LIB_NCMPI   = ${PNETCDF_LIB}
LIB_MPE     =

#LIB_HYPRE   = ${HYPRE_LIB}
LIB_HYPRE   = -I/sw/xt-cle3.1/hypre/2.7.0b/cnl3.1_gnu4.6.1/include -L/sw/xt-cle3.1/hypre/2.7.0b/cnl3.1_gnu4.6.1/lib -lHYPRE

#Specify TEC_PLOT=YES in order to link the tec plot library.
TEC_PLOT=YES
ifeq ($(TEC_PLOT), YES)
CONFIG_LIB = -I${setup_dir}/../source/Simulation/SimulationMain/INavierStokes -L${setup_dir}/../source/Simulation/SimulationMain/INavierStokes -ltecio -lstdc++
endif

#----------------------------------------------------------------------------
# Additional machine-dependent object files
#
#  Add any machine specific files here -- they will be compiled and linked
#  when FLASH is built.
#----------------------------------------------------------------------------

MACHOBJ = 

#----------------------------------------------------------------------------
# Additional commands
#---------------------------------------------------------------------------- 

MV    = mv -f
AR    = ar -r
RM    = rm -f
CD    = cd
RL    = ranlib
ECHO  = echo
AWK   = awk
CAT   = cat
