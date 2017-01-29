./setup INavierStokes/3D/IB_sphere  -3d -auto  -nxb=16 -nyb=16 -nzb=16  -opt -maxblocks=200 -gridinterpolation=native +pm4dev  --without-unit=monitors/Timers -objdir=IB_SPHERE_RE100 PfftSolver=HomBcTrigSolver -site=splash.seas.gwu.edu -noclobber Bittree=1 

