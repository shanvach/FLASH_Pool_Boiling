./setup INavierStokes/3D/SM_OscPlateDEBUG \
  -3d -auto  \
  -nxb=16 -nyb=16 -nzb=16  \
  -debug \
  -maxblocks=200 \
  -gridinterpolation=native \
  +pm4dev  \
  --without-unit=monitors/Timers \
  -objdir=IB_OSCPLATE_RE100_DEBUG \
  PfftSolver=HomBcTrigSolver \
  -site=deepthought.umd.edu_tfitz \
  -noclobber Bittree=1 

