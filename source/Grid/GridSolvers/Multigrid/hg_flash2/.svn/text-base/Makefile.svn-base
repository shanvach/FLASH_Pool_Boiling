# Makefile for Huang-Greengard-style Poisson solver

poisson += hg_level_add.o hg_prolong_bndries.o hg_restrict.o \
           hg_set_ext_bndries.o hg_solve_level.o hg_solve.o \
           poisson_hg_solve_block.o hg_residual.o hg_level_zero.o \
           mg_bndry.o mg_common.o mg_guardcell.o hg_init.o \
           mg_norm.o mg_restore_nodetypes.o mg_set_ext_bndry.o \
           hg_level_smultiply.o hg_init_src.o poisson_image_mass.o \
           poisson_image_boundary.o hg_level_sadd.o

poisson.o hg_set_ext_bndries.o mg_bndry.o mg_guardcell.o \
hg_init.o mg_norm.o mg_restore_nodetypes.o mg_set_ext_bndry.o \
hg_prolong_bndries.o \
hg_residual.o hg_solve.o hg_solve_level.o hg_level_zero.o \
hg_level_add.o hg_level_smultiply.o hg_init_src.o hg_level_sadd.o : mg_common.o
