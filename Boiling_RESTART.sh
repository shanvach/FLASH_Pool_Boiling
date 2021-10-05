./setup INavierStokes/zoso/Boiling_RESTART -2d -auto -nxb=20 -nyb=20 -opt -maxblocks=300 -gridinterpolation=native +pm4dev  -objdir=Boiling_RESTART -site=$1 Bittree=0

git checkout bin/*

git checkout lib/*
