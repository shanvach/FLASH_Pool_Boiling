#	Makefile for multipole Poisson solver

#	Additional files

poisson += init_mpole.o find_center_of_mass.o mpole_common.o \
           compute_mpole_moments.o compute_mpole_potential.o \
           zone_moments.o zone_potential.o

init_mpole.o              : mpole_common.o dBase.o runtime_parameters.o \
                            physical_constants.o logfile.o
find_center_of_mass.o     : mpole_common.o dBase.o
compute_mpole_potential.o : mpole_common.o dBase.o
compute_mpole_moments.o   : mpole_common.o dBase.o
zone_potential.o	  : mpole_common.o dBase.o
zone_moments.o            : mpole_common.o dBase.o

