./setup INavierStokes/zoso/$1 -3d -auto -nxb=16 -nyb=16 -nzb=16 -opt -maxblocks=80 -gridinterpolation=native +pm4dev  -objdir=PoisTest -site=$2 Bittree=0

git checkout bin/*

git checkout lib/*
