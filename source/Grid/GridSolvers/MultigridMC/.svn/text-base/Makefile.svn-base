Grid += $(multigrid)

multigrid = mg_common.o gr_mgInit.o gr_mgInitSlv.o gr_mgInitSrc.o gr_mgNorm.o gr_mgBndry.o \
            gr_mgGuardcell.o gr_mgCycle.o gr_mgZero.o gr_mgCopy.o gr_mgCorrect.o \
            amr_1blk_cc_prol_work_mg.o gr_mgProlong.o gr_mgRestrict.o mg_restore_nodetypes.o  multigrid.o

#mg_correct.o mg_restrict.o mg_cycle.o mg_init.o \
#	      mg_init_src.o \
#             mg_norm.o mg_precond.o mg_prolong.o mg_zero.o multigrid.o \
#	      mg_bndry.o mg_copy.o mg_write.o mg_guardcell.o mg_set_ext_bndry.o \
#             mg_restore_nodetypes.o amr_1blk_cc_prol_work_mg.o

#       New files
#mg_common.o    : dBase.o
#mg_correct.o   : mg_common.o dBase.o
#mg_cycle.o     : mg_common.o dBase.o runtime_parameters.o
#mg_init.o      : mg_common.o dBase.o runtime_parameters.o
#mg_init_src.o  : mg_common.o dBase.o physical_constants.o runtime_parameters.o
#mg_norm.o      : mg_common.o dBase.o
#mg_precond.o   : mg_common.o dBase.o
#mg_prolong.o   : mg_common.o dBase.o
#mg_restrict.o  : mg_common.o dBase.o
#mg_zero.o      : mg_common.o dBase.o
#multigrid.o    : mg_common.o dBase.o runtime_parameters.o logfile.o
#mg_bndry.o     : mg_common.o dBase.o
#mg_copy.o      : mg_common.o dBase.o
#mg_write.o     : mg_common.o dBase.o
#mg_restore_nodetypes.o : mg_common.o dBase.o
#mg_guardcell.o     : mg_common.o dBase.o
#mg_set_ext_bndry.o : mg_common.o dBase.o
