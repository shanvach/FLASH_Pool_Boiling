#-------------------------------------------------------------------
# FLASH makefile definitions for the NERSC IBM SP (seaborg): 64-bit
#
# NOTE: bassi makes a lot of use of Modules, so it is 
# important to make sure that the following modules are
# loaded in order to compile FLASH:
#
#  python
#  
#  hdf5_par_64/1.6.4(default)
#  hdf5_64/1.6.4  (default for serial)
#
# Optional modules:
#
#  GNU
#  MASS
#
# These modules set environment variables pointing to the
# proper library locations. Only python is essential.
#-------------------------------------------------------------------

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------


#HDF5_PATH = /usr/common/usg/hdf5_64/1.4.5-post2/parallel
HDF5_PATH = /usr/common/usg/hdf5_64/1.6.4/parallel
#HDF5_PATH = /usr/common/usg/hdf5_64/1.6.4/serial
ZLIB_PATH = /usr/common/usg/zlib_64/1.2.1/lib

#HDF5 was built with the SZIP library enabled.  Need to link to it
SZIP_PATH = /usr/common/usg/szip/2.0/lib

PNG_PATH  = /usr/common/usg/gnu

HPM_PATH  = /usr/common/usg/hpmtoolkit/2.5.2

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

FCOMP   = mpxlf90_r
CCOMP   = mpcc_r
CPPCOMP = mpCC_r
LINK    = mpxlf90_r

#----------------------------------------------------------------------------
# Compilation flags
#
#  Three sets of compilation/linking flags are defined: one for optimized
#  code, one for testing, and one for debugging.  The default is to use the 
#  _OPT version.  Specifying -debug to setup will pick the _DEBUG version,
#  these should enable bounds checking.  Specifying -test is used for 
#  flash_test, and is set for quick code generation, and (sometimes) 
#  profiling.  The Makefile generated by setup will assign the generic token 
#  (ex. FFLAGS) to the proper set of flags (ex. FFLAGS_OPT).
#----------------------------------------------------------------------------

FFLAGS_OPT   = -O2 -qintsize=4 -qrealsize=8 -qfixed -qnosave -q64 -c \
               -qinline -qnoipa -qmaxmem=16384 -qxlf90=autodealloc \
               -qsuffix=cpp=F -qarch=auto -qtune=auto -qcache=auto 

FFLAGS_TEST  = -O2 -qintsize=4 -qrealsize=8 -qfixed -qnosave -q64 -c \
               -qsuffix=cpp=F -qarch=auto -qtune=auto -qcache=auto 

FFLAGS_DEBUG = -g -qintsize=4 -qrealsize=8 -qfixed -qnosave -q64 -c \
               -qarch=auto 


F90FLAGS     = -qsuffix=f=F90:cpp=F90 -qfree=f90
f90FLAGS     = -qsuffix=f=f90:cpp=F90 -qfree=f90



# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5  = -I${HDF5_PATH}/include

CFLAGS_OPT   = -O2 -DIBM -DNOUNDERSCORE -q64 -c \
               -qarch=auto -qtune=auto -qcache=auto
CFLAGS_TEST  = -O2 -DIBM -DNOUNDERSCORE -q64 -c \
               -qarch=auto -qtune=auto -qcache=auto
CFLAGS_DEBUG = -g  -DIBM -DNOUNDERSCORE -q64 -c \
               -qarch=auto

MDEFS = -WF,

  .SUFFIXES: .o .c .f .F .h .fh .F90 .f90

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------



LFLAGS_OPT   = -q64 -o
LFLAGS_TEST  = -q64 -o
LFLAGS_DEBUG = -q64 -g -o

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

LIB_MPI =

  LIB_HDF5  = -L${HDF5_PATH}/lib -lhdf5 \
  -L${ZLIB_PATH} -lz \
  -L${PNG_PATH}/lib -lpng \
  -L${SZIP_PATH} -lsz -lgpfs

LIB_MATH  = -lessl

LIB_OPT   = ${MASS}
LIB_DEBUG =
LIB_TEST  =

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

MV = mv -f
AR = ar -r
RM = rm -f
CD = cd
RL = ranlib
ECHO = echo
