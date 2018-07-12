./setup INavierStokes/zoso/Nucleate_Boiling_3D_tests  -3d -auto  -nxb=20 -nyb=20 -nzb=20 -opt  -maxblocks=16 -gridinterpolation=native +pm4dev  -objdir=Nucleate_Boiling_3D -parfile=$1 -site=$2 Bittree=0

git checkout bin/*

git checkout lib/*
