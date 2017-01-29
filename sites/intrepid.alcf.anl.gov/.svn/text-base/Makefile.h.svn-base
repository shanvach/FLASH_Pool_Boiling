#-------------------------------------------------------------------
# FLASH makefile definitions for ALCF (argonne) production BG/P (intrepid)
#  started from seaborg (NERSC) makefile
#
#
#-------------------------------------------------------------------

MASS = -lmass

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

HDF5_PATH = /soft/apps/hdf5-1.8.0

ZLIB_PATH =

HPM_PATH =
PNG_PATH =

NCMPI_PATH = /bgusr/robl/soft/pnetcdf-20071019
MPI_PATH =

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

FCOMP   = mpixlf90_r
CCOMP   = mpixlc_r
CPPCOMP = mpixlcxx_r
LINK    = mpixlf90_r

#This uses an old version of hpctoolkit which is known to be buggy.
#See https://mailman.rice.edu/pipermail/hpctoolkit-forum/2011-October/000301.html
#soft add +hpctoolkit
#LINK    = hpclink mpixlf90_r

#This is a much newer version of hpctoolkit
#LINK = /home/projects/hpctoolkit/pkgs/hpctoolkit/bin/hpclink mpixlf90_r

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

#See special BG/P OpenMP compilation at the bottom of this file.
OPENMP = -qsmp=omp:noauto
FFLAGS_OPT   = -g -O3 -qnohot -qintsize=4 -qrealsize=8 -qfixed -qnosave -c \
               -qarch=450 -qtune=auto -qcache=auto

FFLAGS_TEST  = -g -O3 -qnohot -qintsize=4 -qrealsize=8 -qfixed -qnosave -c \
               -qarch=450 -qtune=auto -qcache=auto

#Use -O2 because -O0 generates really slow code.
FFLAGS_DEBUG = -g -O2 -qintsize=4 -qrealsize=8 -qfixed -qnosave -c \
               -qarch=450 -qtune=auto -qcache=auto \
               -qfloat=rngchk -qcheck -qflttrap=enable:invalid:nanq:overflow:zerodivide

F90FLAGS     = -I${HDF5_PATH}/include -DH5_USE_16_API -qsuffix=f=F90:cpp=F90 -qfree=f90
f90FLAGS     = ${F90FLAGS}

# if we are using HDF5, we need to specify the path to the include files

CFLAGS_OPT   = -g -O3 -qnohot -DIBM -c -D_FILE_OFFSET_BITS=64 \
               -qarch=450 -qtune=auto -qcache=auto
CFLAGS_TEST  = -g -O3 -qnohot -DIBM -c -D_FILE_OFFSET_BITS=64 \
               -qarch=450 -qtune=auto -qcache=auto
CFLAGS_DEBUG = -g -O2 -DIBM -c -D_FILE_OFFSET_BITS=64 \
               -qarch=450 -qtune=auto -qcache=auto -qdbxextra \
               -qfloat=rngchk -qcheck=all -qflttrap=enable:invalid:nanq:overflow:zerodivide

#Add include path to CFLAGS_MPI for the memory sampler in ut_sysMemBGKernel.c
CFLAGS_MPI = -I/bgsys/drivers/ppcfloor/arch/include
CFLAGS_HDF5  = -I${HDF5_PATH}/include -DH5_USE_16_API
CFLAGS_NCMPI = -I$(NCMPI_PATH)/include

MDEFS = -WF,

.SUFFIXES: .o .c .f .F .h .fh .F90 .f90

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT,
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -g -O3 -qnohot -o
LFLAGS_TEST  = -g -pg -O3 -qnohot -o
LFLAGS_DEBUG = -g -O2 -o

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
#LIB_MPI   = -L/home/morozov/lib -lmpihpm
#LIB_MPI   = -L /soft/apps/current/ibm-hpct/lib -lmpitrace -llicense

#No need to distinguish between mpihpm libraries on BG/P.
LIB_MPIHPM = -L/soft/apps/current/hpm/lib -lmpihpm
LIB_MPIHPM_SMP = ${LIB_MPIHPM}

LIB_HDF5  = -L${HDF5_PATH}/lib -lhdf5_fortran -lhdf5
LIB_NCMPI = -L$(NCMPI_PATH)/lib -lpnetcdf
LIB_LAPACK = -L/soft/apps/LAPACK -llapack_bgp -L/soft/apps/BLAS -lblas_bgp

LIB_MATH  =

LIB_OPT   = ${MASS}
LIB_DEBUG =
LIB_TEST  = ${MASS}

LIB_STDCXX = -L$(shell dirname $(shell dirname $(XLC_USR_CONFIG)))/lib -libmc++ -lstdc++

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

# This next section only applies to compiling FLASH, not some library.
# Anything that mentions a specific file should be within this ifeq block
ifeq ($(FLASHBINARY),true)

# add_block_to_tree subroutine in local_tree.F90 is called with the same argument twice.
# The arguments alias each other so we have to inform the compiler not to use unsafe
# optimisations (-qalias).  Required at -O4 optimisation level for correctness.
local_tree_module.mod local_tree.mod local_tree.o : local_tree.F90
	@$(ECHO) Specially compile $< with -qalias=nostd flag
	${FCOMP} ${FFLAGS} -qalias=nostd ${F90FLAGS} ${FDEFINES} $<

endif
#End of FLASH binary if statement.
