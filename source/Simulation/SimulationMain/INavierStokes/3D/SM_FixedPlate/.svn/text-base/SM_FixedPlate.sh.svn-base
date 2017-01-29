#./setup SolidMechanics/Plate -3d -auto -nxb=4 -nyb=4 -nzb=4 -debug -maxblocks=1 +pm4dev -objdir=SolidMechanics_Plate -site=crash2.umd.edu -noclobber Bittree=1 -gridinterpolation=native

./setup INavierStokes/3D/SM_FixedPlate \
        -3d \
        -auto \
        -debug \
        -nxb=16 \
        -nyb=16 \
        -nzb=16 \
        -maxblocks=300 \
        -gridinterpolation=native \
        +pm4dev \
        --without-unit=monitors/Timers \
        -objdir=SM_FixedPlate \
         PfftSolver=HomBcTrigSolver \
        -site=deepthought.umd.edu_tfitz \
        -noclobber \
        Bittree=1 
