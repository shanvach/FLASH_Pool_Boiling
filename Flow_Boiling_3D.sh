./setup INavierStokes/zoso/Flow_Boiling_3D  -3d -auto  -nxb=20 -nyb=20 -nzb=20 -opt  -maxblocks=300 -gridinterpolation=native +pm4dev  -objdir=Flow_Boiling_3D -site=$1 Bittree=0

git checkout bin/*

git checkout lib/*
