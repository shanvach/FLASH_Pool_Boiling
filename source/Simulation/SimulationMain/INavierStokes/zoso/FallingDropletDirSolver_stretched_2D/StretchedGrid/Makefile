
#Location of HDF5 binaries (with include/ and lib/ underneath)
HDF5_PATH = /usr
MPI_PATH = /usr/lib64/mpich

# Compiler
FC   = ${MPI_PATH}/bin/mpif90

# ------ No machine-specific paths/variables after this  -----

LIB_HDF5  = -L${HDF5_PATH}/lib -lhdf5_fortran -lhdf5 -lz
LIB_MPI = -L${MPI_PATH}/lib

FSOURCE = gridGen
OBJECTS = gridGen.o

FLAGS   = -c

grid:  gridGen.F90
	$(FC) -o $@ gridGen.F90 $(LIB_HDF5) $(LIB_MPI)

clean:
	rm -f $(FSOURCE) $(OBJECTS) *.h5

.F90.o: $(FSOURCE)
	$(FC) $(FLAGS) $?

.SUFFIXES:.o.c.F90
