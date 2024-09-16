./setup INavierStokes/zoso/Pool_Boiling_2D  -2d -auto  -nxb=16 -nyb=16  -opt  -maxblocks=2000 -gridinterpolation=native +pm4dev  -objdir=Pool_Boiling_2D -site=$1 Bittree=0

git checkout bin/*

git checkout lib/*
