#       Makefile.h file for cetus.asci.uchicago.edu
#
#	FLASH makefile definitions for Linux (Lahey compiler)


#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------


HDF5_PATH   = /opt/hdf5/lahey/1.8.7
MPI_PATH    = /opt/mpich2/lahey/1.4.1p1


ZLIB_PATH  =

NCMPI_PATH = /opt/netcdf/lahey/1.2.0
MPE_PATH   =
HYPRE_PATH = /opt/hypre/lahey/2.8.0b
ifeq ("$(USEOPENMP)", "1")
HYPRE_PATH=/opt/hypre/lahey/2.8.0b_omp
endif

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

FCOMP      = $(MPI_PATH)/bin/mpif90 
CCOMP      = $(MPI_PATH)/bin/mpicc
CPPCOMP    = $(MPI_PATH)/bin/mpicxx
LINK       = $(MPI_PATH)/bin/mpif90


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

#I've noticed on our VM test platform that the number of OpenMP threads
#cannot be greater than the number of processors on the machine, i.e.
#setting OMP_NUM_THREADS to 4 or higher results in only 3 OpenMP threads.

#The -Kpureomp option forces strict adherence to OpenMP directives.
#-Knopureomp allows the compiler to optimize OpenMP code due to generous
#interpretation of OpenMP directives
OPENMP_FORTRAN = --openmp -Kpureomp
OPENMP_C = -fopenmp
OPENMP_LINK = --openmp -Kpureomp

FFLAGS_OPT   = -c --o2 -CcdRR8
FFLAGS_DEBUG = -c -g --trace --trap --chk[aes] -CcdRR8

#The FFLAGS_OPT optimization is switched off because the MPI wrapper
#script appends a '-g' to the flags - therefore in FFLAGS_TEST I turn off
#the debugging information.  I also ensure that we get consistent
#answers between code compiled with and without OpenMP by using the
#'--ap' flag.  This is based on a recent experience with Intel compiler 
#where the floating point transformations differed between code
#compiled (and optimized) with and without OpenMP support.
# -[n]g                    generate debugging information
#--[n]ap                   ensure consistent arithmetic precision
FFLAGS_TEST  = ${FFLAGS_OPT} -ng --ap
#FFLAGS_TEST  = -c -g --trace --trap --chk[aesx] --chkglobal -CcdRR8

FFLAGS_HYPRE = -I${HYPRE_PATH}/include
CFLAGS_HYPRE = -I${HYPRE_PATH}/include


F90FLAGS     =

#CFLAGS       = -c -O3 -tpp7 -march=pentium4 -mcpu=pentium4 -ip -unroll \
#               -D_LARGEFILE64_SOURCE

CFLAGS_OPT =   -I$(MPI_PATH)/include -c -O2 -D_LARGEFILE64_SOURCE -g
CFLAGS_DEBUG = -I$(MPI_PATH)/include -c -g 
CFLAGS_TEST  = -I$(MPI_PATH)/include -c -g
CFLAGS_HDF5  = -I$(HDF5_PATH)/include -DH5_USE_16_API
CFLAGS_NCMPI = -I$(NCMPI_PATH)/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT = -o
LFLAGS_TEST = -ng -o
LFLAGS_DEBUG = -o  

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

LIB_MPI     =
LIB_HDF5    = -L$(HDF5_PATH)/lib -lhdf5 -lz

LIB_PAPI    = #$(PAPI_PATH)/lib/libpapi.a $(PAPI_PATH)/lib/_fixunssfdi.o
LIB_PNG     = #-lpng -lz

LIB_OPT     =
LIB_DEBUG   =
LIB_TEST    =

LIB_NCMPI   = -L$(NCMPI_PATH)/lib -lpnetcdf

LIB_MPE     =
LIB_HYPRE = -L${HYPRE_PATH}/lib -lHYPRE
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


#----------------------------------------------------------------------------
# Specific compilation.
#  fftpack mixes reals and integers in 
#  procedure calls.  We tell the compiler to compile this file
#  with more forgiving runtime checking.
#---------------------------------------------------------------------------- 
ifeq ($(FLASHBINARY),true)

FFLAGS_WO_CHKA = $(patsubst --chk[a%],--chk[%],$(FFLAGS))
FFLAGS_WO_CHKG = $(filter-out --chkglobal,$(FFLAGS_WO_CHKA))
# Compile the following file with the same flags as others (which depend on
# whether -opt, -debug, or -test was in effect for setup), except that
# --chk[aes] is replaced by --chk[es]. This allows compilation of this file
# with debugging as used on cetus with the Lahey compiler.
#
# Without this workaround, the code fails at runtime with messages like this:
#  The type of argument 3 is inconsistent (actual argument wsave: r*8, dummy argument ifac: i*4).
# in rffti1.
#
fftpack.o : %.o : %.f90
	$(FCOMP) $(FFLAGS_WO_CHKG) $(F90FLAGS) $(FDEFINES) $<

FFLAGS_WO_CHKU0 = $(patsubst --chk[%ux],--chk[%x],$(FFLAGS))
FFLAGS_WO_CHKU = $(patsubst --chk[%u],--chk[%],$(FFLAGS_WO_CHKU0))
gr_bcGetRegion.o : %.o : %.F90 physicaldata.o Driver_interface.o constants.h Flash.h
	$(FCOMP) $(FFLAGS_WO_CHKU) $(F90FLAGS) $(FDEFINES) $<
endif

