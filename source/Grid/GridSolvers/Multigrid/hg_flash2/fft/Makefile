# Additional files and dependencies for FFT-based coarse-grid solver

poisson += pois_solve_1d.o pois_solve_2d.o pois_solve_3d.o \
           fftsg.o fftsg2d.o fftsg3d.o

poisson_hg_solve_block.o : pois_solve_1d.o pois_solve_2d.o pois_solve_3d.o \
                           dBase.o
