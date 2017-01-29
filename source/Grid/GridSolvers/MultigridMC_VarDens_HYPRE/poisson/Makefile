#	Makefile for multigrid Poisson solver

Grid += Grid_solvePoisson.o poisson_mg_residual.o poisson_mg_residualMG.o \
  poisson_mg_relax_ZEBRA.o \
  poisson_mg_relax_RBGS.o \
  poisson_mg_relax.o \
  poisson_mg_solve.o poisson_mg_relax_RBGSMG.o \
  poisson_mg_relax_HYPRE.o gr_mgCycleHYPRE.o \
  gr_mgCycleHYPRE_PC.o \
  multigridHYPRE.o multigridHYPRE_PC.o


#          poisson_mg_relax.o poisson_mg_residual.o \
#          poisson_mg_solve.o 


#       New files

#poisson_mg_relax.o     : mg_common.o \
#	                  dBase.o     \
#		          physical_constants.o \
#                         runtime_parameters.o \
#                         logfile.o
#poisson_mg_residual.o  : mg_common.o dBase.o
#poisson_mg_solve.o     : mg_common.o dBase.o runtime_parameters.o

#poisson_image_mass.o : dBase.o

#	Additional dependencies

#poisson.o      : multigrid.o mg_common.o dBase.o logfile.o

