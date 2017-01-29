./setup INavierStokes/3D/RBC_SHEAR  -3d -auto  -nxb=16 -nyb=16 -nzb=16  -opt -maxblocks=100 -gridinterpolation=native +pm4dev  --without-unit=monitors/Timers -objdir=SM_SHEAR_RE100 PfftSolver=HomBcTrigSolver -site=SEAS12982.gwu.edu -noclobber Bittree=1 

