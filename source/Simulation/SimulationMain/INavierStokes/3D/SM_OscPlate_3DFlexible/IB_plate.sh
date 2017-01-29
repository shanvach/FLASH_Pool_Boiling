./setup INavierStokes/3D/SM_OscPlate_3DFlexible \
  -3d -auto  \
  -nxb=16 -nyb=16 -nzb=16  \
  -opt \
  -maxblocks=200 \
  -gridinterpolation=native \
  +pm4dev  \
  --without-unit=monitors/Timers \
  -objdir=IB_OSCPLATE_RE100_3DFlex \
  PfftSolver=HomBcTrigSolver \
  -site=deepthought.umd.edu_tfitz \
  -noclobber Bittree=1 

