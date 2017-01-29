./setup INavierStokes/2D/FallCylinder/  -2d -auto -opt +ug -nxb=2640 -nyb=66 -maxblocks=1 -gridinterpolation=native PfftSolver=HomBcTrigSolver  -objdir=MOVECYL_2D_48D_PRESCRIBED -site=colonialone.gwu.edu


#./setup INavierStokes/2D/TaylorGreenVortex  -2d -auto +ug -nxb=64 -nyb=32 -debug -objdir=TAYLOR_GREEN_UG_GEN -site=SEAS10927.gwu.edu -noclobber

#./setup INavierStokes/2D/TaylorGreenVortex  -2d -auto  -nxb=128 -nyb=128 -opt  -maxblocks=1  -parfile=/Users/mvanella/Documents/ADAPTIVE/FLASH/test/INS_IB/source/Simulation/SimulationMain/INavierStokes/2D/TaylorGreenVortex/flash.par +ug  -objdir=TAYLOR_GREEN_IB_CYL PfftSolver=HomBcTrigSolver -site=SEAS10927.gwu.edu -noclobber 

# -gridinterpolation=native
# +pm4dev
# PfftSolver=HomBcTrigSolver
