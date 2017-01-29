# FLASH makefile definitions for ix86 Linux (Portland Group compiler)
#
# Red Hat 9
# Portland Group Fortran 5.2-4
# hdf5-1.6.2 (gcc-3.2.2)
# mpich-1.2.6 (gcc 3.2.2)

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

HDF5_PATH  = /share/home/01173/jzuhone/opt
ZLIB_PATH  =

MPE_PATH   =

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

FCOMP      = mpif90
CCOMP      = mpicc
CPPCOMP    = mpicxx
LINK       = mpif90

# pre-processor flag

PP         = -D

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

FFLAGS_OPT   = -c -r8 -i4 -O3
FFLAGS_DEBUG = -c -r8 -i4 -g -Ktrap=divz -Mchkfpstk -Mchkptr -Mchkstk -Mdclchk -Mbounds
FFLAGS_TEST  = -c -r8 -i4 -Mprof=lines

CFLAGS_OPT   = -c -O2
CFLAGS_DEBUG = -c -g
CFLAGS_TEST  = -c 

# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5  = -I${HDF5_PATH}/include
CFLAGS_NCMPI = -I${TACC_NETCDF_INC}
CFLAGS_GSL   = -I${TACC_GSL_INC} -I${TACC_GSL_INC}/gsl

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -o
LFLAGS_DEBUG = -Bstatic -Wl,-noinhibit-exec -g -o
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


LIB_HDF5    = -L${HDF5_PATH}/lib -lhdf5 -lz
              
LIB_MPI     = 
LIB_PNG     = -lpng -lz

LIB_OPT     = 
LIB_DEBUG   =
LIB_TEST    =

LIB_GSL     = -L${TACC_GSL_LIB} -lgsl
LIB_NCMPI   = -L${TACC_NETCDF_LIB} -lnetcdf
LIB_MPE     =

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


