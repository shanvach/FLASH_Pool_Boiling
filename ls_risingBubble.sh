./setup INavierStokes/2D/LS_rising_bubble  -2d -auto  -nxb=20 -nyb=20 -opt  -maxblocks=10 -gridinterpolation=native +pm4dev  -objdir=LS_RISING_BUBBLE -site=$1 Bittree=0

git checkout bin/*

git checkout lib/*
