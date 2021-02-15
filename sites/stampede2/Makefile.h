# FLASH makefile definitions for x86-64 Linux
#
#----------------------------------------------------------------------------
# Set the HDF5/MPI library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

HYPRE_PATH   =
SUPERLU_PATH = 

#----------------------------------------------------------------------------
# Compiler and linker commands
#----------------------------------------------------------------------------

FCOMP   = ${TACC_IMPI_BIN}/mpif90
CCOMP   = ${TACC_IMPI_BIN}/mpicc
CPPCOMP = ${TACC_IMPI_BIN}/mpicxx
LINK    = ${TACC_IMPI_BIN}/mpif90

# pre-processor flag
PP      = -D

#----------------------------------------------------------------------------
# Compilation flags
#----------------------------------------------------------------------------

FFLAGS_OPT = -c -r8 -i4 -real_size 64 ${TACC_VEC_FLAGS} -align array32byte -O3
FFLAGS_DEBUG = -c -r8 -i4 -real_size 64 -O0 -check bounds -check format \
               -check output_conversion -warn error -check uninit -traceback -fp-stack-check -fpe0 -check pointers
FFLAGS_TEST  = ${FFLAGS_OPT} -fp-model precise

F90FLAGS = -DH5_USE_16_API -D_LARGEFILE64_SOURCE -D_FORTIFY_SOURCE=2 

CFLAGS_OPT = -c ${TACC_VEC_FLAGS} -O3 -D_LARGEFILE64_SOURCE
CFLAGS_DEBUG = -c -O0 -g -traceback -debug all -debug extended \
               -D_LARGEFILE64_SOURCE -ftrapuv -fp-stack-check
CFLAGS_TEST = ${CFLAGS_OPT} -fp-model precise

CFLAGS_HDF5 = -I${TACC_HDF5_INC} -DH5_USE_16_API
FFLAGS_HDF5 = -I${TACC_HDF5_INC} -DH5_USE_16_API

CFLAGS_HYPRE = -I${HYPRE_PATH}/include
FFLAGS_HYPRE = -I${HYPRE_PATH}/include

CFLAGS_SUPERLU = -I${SUPERLU_PATH}/include
FFLAGS_SUPERLU = -I${SUPERLU_PATH}/include

#----------------------------------------------------------------------------
# Linker flags
#----------------------------------------------------------------------------

LFLAGS_OPT   = -O3 -o
LFLAGS_DEBUG = -o
LFLAGS_TEST  = -o

#----------------------------------------------------------------------------
# Library specific linking
#----------------------------------------------------------------------------

LIB_OPT   = 
LIB_DEBUG = 
LIB_TEST  =

LIB_HDF5    = -Wl,-rpath,${TACC_HDF5_LIB} -L${TACC_HDF5_LIB} -lhdf5_fortran -lhdf5 -lz
LIB_MPI     = 
LIB_BLAS    = -mkl
LIB_HYPRE   = -L${HYPRE_PATH}/lib -lHYPRE  
LIB_SUPERLU = -L${SUPERLU_PATH}/lib -lsuperlu_4.3
LIB_STDCXX  = -L/lib64 -lstdc++

#----------------------------------------------------------------------------
# Additional machine-dependent object files
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

