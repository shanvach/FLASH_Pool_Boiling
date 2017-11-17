# FLASH makefile definitions for the Cray XT5 System, Kraken at NICS
#
# The XT5 makes use of modules.
# HDF5: module load hdf5-parallel
# PNETCDF: module load p-netcdf
#
# You should load all required modules before compiling!
# See http://www.nics.tennessee.edu/computing-resources/kraken/compiling
# In this Makefile.h we select the compilation flags based on the
# programming environment.  If you need HDF5 then you should load HDF5
# before compiling because we need the compiler wrapper scripts to
# include the HDF5 header files and also link against the HDF5 library.
#
#----------------------------------------------------------------------------
# Set the HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

export cur-dir := $(shell pwd)

# Set the location of top directory
export setup_dir = $(cur-dir)

#PE_ENV=GNU
PE_ENV=INTEL

MPI_PATH   =
PAPI_PATH  =
PAPI_FLAGS =
NCMPI_PATH =
HYPRE_PATH = /p/home/balaras/software/hypre_library
SUPERLU_PATH = 
LAPACK_PATH =
BLAS_PATH = 
#BLAS_PATH = ${BLAS_DIR}
#HDF5_PATH = ${HDF5DIR}
HDF5_PATH = /p/home/balaras/software/hdf5_library
SUPERLUDIST_PATH = 

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

ifeq ($(findstring $(PE_ENV),PGI GNU CRAY INTEL),)
$(error Your environment "$(PE_ENV)" is invalid.  It must be "PGI", "GNU" or "CRAY")
else
$(warning You are using the "$(PE_ENV)" environment)
endif

ifdef CRAY_HDF5_DIR
$(warning You are using the HDF5 installation at "$(CRAY_HDF5_DIR)")
else
$(warning Generally FLASH needs HDF5.  You can load it using "module load hdf5-parallel")
endif


#FCOMP   =  gfortran
#CCOMP   =  gcc
#CPPCOMP =  g++
#LINK    =  gfortan

#FCOMP   =  ifort
#CCOMP   =  icc
#CPPCOMP =  icpc
#LINK    =  ifort

#FCOMP   =  mpif90
#CCOMP   =  mpicc
#CPPCOMP =  mpicxx
#LINK    =  mpif90

FCOMP   = ftn
CCOMP   = cc
CPPCOMP = CC
LINK    = ftn

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

# PGI Compiler:
PGI_FFLAGS_OPT   = -c -r8 -i4 -fastsse -Minfo -Mneginfo -tp amd64
PGI_FFLAGS_DEBUG = -c -r8 -i4 -g -Ktrap=divz -Mchkfpstk -Mchkptr -Mchkstk -Mdclchk -Mbounds
PGI_FFLAGS_TEST  = -c -r8 -i4 -Mprof=lines

# GNU Compiler:
GNU_FFLAGS_OPT =  -c -O2 -fdefault-real-8 -fdefault-double-8 \
-ffree-line-length-none -Wuninitialized
GNU_FFLAGS_DEBUG = -ggdb -c -fdefault-real-8 -fdefault-double-8 \
-ffree-line-length-none -pedantic -Wall -Wextra -Waliasing \
-Wsurprising -Wconversion -Wunderflow \
-ffpe-trap=invalid,zero,overflow -fbounds-check \
-fbacktrace -fdump-core -finit-real=nan \
-finit-integer=-999999 -fimplicit-none -Wstack-protector -g -dynamic
GNU_FFLAGS_TEST =  -c -g -dynamic -fdefault-real-8 -fdefault-double-8 \
-ffree-line-length-none

# CRAY Compiler:
CRAY_FFLAGS_OPT   = -c -O2 -s real64 -s integer32
CRAY_FFLAGS_DEBUG = -c -G 0 -s real64 -s integer32
CRAY_FFLAGS_TEST  = -c -s real64 -s integer 32

# INTEL compiler:
INTEL_FFLAGS_OPT   = -c -O3 -r8 -i4 -real_size 64
#-check uninit
#-optTauSelectFile="/nics/b/home/mvanella/ADAPTATIVE/FLASH4_INS/sites/kraken.nics.tennessee.edu/select.tau"
INTEL_FFLAGS_DEBUG = -c -g -r8 -i4 -O0 -check bounds -check format \
-check output_conversion -warn error -real_size 64 -check uninit \
-traceback -fp-stack-check  -fpe0 -check pointers
INTEL_FFLAGS_TEST  = -c -g -r8 -i4


#Now we use the correct compiler flags:
FFLAGS_OPT = $($(PE_ENV)_FFLAGS_OPT)
FFLAGS_DEBUG = $($(PE_ENV)_FFLAGS_DEBUG)
FFLAGS_TEST = $($(PE_ENV)_FFLAGS_TEST)

CFLAGS_OPT   = -c -O3
CFLAGS_DEBUG = -c -O1 -g -dynamic
CFLAGS_TEST  = -c

#No path required because we are using compiler wrapper scripts.
#FFLAGS_MPI   = -I${CRAY_MPICH2_DIR}/include
#CFLAGS_MPI   = -I${CRAY_MPICH2_DIR}/include

#FFLAGS_HDF5 = -I${HDF5INCLUDE} -DH5_USE_16_API  -I${HDF5DIR} -DH5_USE_16_API
FFLAGS_HDF5 = -I${HDF5_PATH}/include -DH5_USE_16_API
#CFLAGS_HDF5 = -I${HDF5INCLUDE} -DH5_USE_16_API  -I${HDF5DIR} -DH5_USE_16_API
CFLAGS_HDF5 = -I${HDF5_PATH}/include -DH5_USE_16_API
CFLAGS_NCMPI =

FFLAGS_PAPI  =

CFLAGS_SUPERLU = -I${SUPERLU_PATH}/SRC
FFLAGS_SUPERLU = -I${SUPERLU_PATH}/SRC

CFLAGS_PARMETIS = -I${PARMETIS_PATH}/include
FFLAGS_PARMETIS = -I${PARMETIS_PATH}/include

CFLAGS_HYPRE = -I${HYPRE_PATH}/include
FFLAGS_HYPRE = -I${HYPRE_PATH}/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT,
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

ifeq ($(PE_ENV), INTEL)

LFLAGS_OPT   = -diag-disable 10120 -O3 -o
LFLAGS_DEBUG = -diag-disable 10120 -o
LFLAGS_TEST  = -diag-disable 10120 -O3 -o

else

LFLAGS_OPT   = -O3 -o
LFLAGS_DEBUG = -o
LFLAGS_TEST  = -O3 -o

endif

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

#LIB_HDF5    = -L${HDF5DIR} -lhdf5_fortran -lhdf5
LIB_HDF5    = -L${HDF5_PATH}/lib -lhdf5_fortran -lhdf5 -lz
LIB_MPI     = 
LIB_PAPI    =
LIB_PNG     =

ifeq ($(PE_ENV), INTEL)
LIB_MATH  = -limf -lm
else
LIB_MATH  = -lm
endif

LIB_OPT     =
LIB_DEBUG   =
LIB_TEST    =

LIB_NCMPI   =
LIB_MPE     =

LIB_BLAS    = -L${BLAS_PATH} -lblas_LINUX
LIB_SUPERLU = -L${SUPERLU_PATH}/lib -lsuperlu_4.3
LIB_SUPERLUDIST  = 
LIB_LAPACK  =
LIB_HYPRE   = -L${HYPRE_PATH}/lib -lHYPRE 

#LIB_STDCXX  = -L/usr/lib -lstdc++
LIB_STDCXX  = -lstdc++

#Specify TEC_PLOT=YES in order to link the tec plot library.
TEC_PLOT=YES
ifeq ($(TEC_PLOT), YES)
# CONFIG_LIB = -I/usr/people/brownmj/Software/tecio -L/usr/people/brownmj/Software/tecio -ltecio -lstdc++
CONFIG_LIB = -I/p/home/balaras/software -L/p/home/balaras/software -ltecio -lstdc++
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
