./setup INavierStokes/zoso/Pool_Boiling_2D -2d -auto -nxb=20 -nyb=20 -opt -maxblocks=300 -gridinterpolation=native +pm4dev  -objdir=Pool_Boiling_2D -site=$1 Bittree=0

git checkout bin/*

git checkout lib/*
