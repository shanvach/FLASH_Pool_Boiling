# FLASH makefile definitions for the compilers on the ramses cluster

#OPENMP = -fopenmp

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

HDF4_PATH=${HDF_HOME}
HDF5_PATH=/home/orban/hdf5-1.8.7_parallel/hdf5

ZLIB_PATH  =

PAPI_PATH  =
PAPI_FLAGS =

NCMPI_PATH =
MPE_PATH   =

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

FCOMP   = mpif90
CCOMP   = mpicc
CPPCOMP = mpiCC
LINK    = mpif90

# pre-processor flag
PP     = -D

#----------------------------------------------------------------------------
# Compilation flags
#
#  Three sets of compilation/linking flags are defined: one for optimized
#  code, one for testing, and one for debugging.  The default is to use the 
#  _OPT version.  Specifying -debug to setup will pick the _DEBUG version,
#  these should enable bounds checking.  Specifying _TEST is used for 
#  flash_test, and is set for quick code generation, and (sometimes) 
#  profiling.  The Makefile generated by setup will assign the generic token 
#  (ex. FFLAGS) to the proper set of flags (ex. FFLAGS_OPT).
#----------------------------------------------------------------------------

FFLAGS_OPT   = -c -O2 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -I/home/orban/hypre-2.7.0b/src/
FFLAGS_DEBUG = -c -g -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -I/home/orban/hypre-2.7.0b/src/
FFLAGS_TEST  = -c -fdefault-real-8 -fdefault-double-8

CFLAGS_OPT   = -O2 -c
CFLAGS_DEBUG = -g -c
CFLAGS_TEST  = -c

CFLAGS_HDF5  = -I$(HDF5_PATH)/include -DH5_USE_16_API
CFLAGS_NCMPI =

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -o
LFLAGS_DEBUG = -g -o
LFLAGS_TEST  = -o


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

#LIB_HDF4  = -L${HDF4_PATH}/lib -lmfhdf -ldf -lz -ljpeg
LIB_HDF5  = -L${HDF5_PATH}/lib -lhdf5 -lz

LIB_HYPRE = -L/home/orban/hypre-2.7.0b/src/lib -lHYPRE

LIB_OPT   = 
LIB_DEBUG = 
LIB_TEST  =

LIB_NCMPI =
LIB_MPE   =
LIB_MPI   =

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

MV   = mv -f
AR   = ar -r
RM   = rm -f
CD   = cd
RL   = ranlib
ECHO = echo



