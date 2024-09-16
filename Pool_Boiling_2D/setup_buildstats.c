/* !!! DO NOT EDIT, FILES WRITTEN BY SETUP SCRIPT !!!
   
!!****f* object/setup_buildstats
!!
!! NAME
!!
!!  setup_buildstats
!!
!!
!! SYNOPSIS
!!
!!  call setup_buildstats(build_date, build_dir, build_machine, setup_call)
!!
!!  call setup_buildstats(character, character, character, character)
!!
!! 
!! DESCRIPTION
!!
!!  Simple subroutine generated at build time that returns the build date
!!  build directory, build machine, c and f flags and the full setup 
!!  command used to assemble the FLASH executable.
!!
!!
!!
!!***
*/

#include "Flash.h"
#include "constants.h"
#include "mangle_names.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

void FTOC(setup_buildstats)(char* build_date, 
		    char* build_dir, 
		    char* build_machine, 
		    char* setup_call, 
		    char* c_flags, 
		    char* f_flags){



     strncpy(build_date, "Wed Apr 24 15:53:34 PDT 2024",80);
     strncpy(build_dir, "/home3/svachhan/FLASH/Pool_Boiling_2D", 80);
     strncpy(build_machine, "Linux pfe20 4.18.0-477.27.1.1toss.t4.x86_64 #1 SMP Tue Sep 19 15:17:56 PDT 2023 ", 80);
     strncpy(setup_call, "/home3/svachhan/FLASH/bin/setup.py INavierStokes/zoso/Pool_Boiling_2D -2d -auto -nxb=20 -nyb=20 -opt -maxblocks=2000 -gridinterpolation=native +pm4dev -objdir=Pool_Boiling_2D -site=pfe.nas.nasa.gov Bittree=0 ",400);
     strncpy(c_flags, "/nasa/sgi/mpt/2.15r20/bin/mpicc -I/home3/svachhan/hdf5-1.8.20_NEW/hdf5/include -DH5_USE_16_API -I/home3/svachhan/hypre-2.11.2_NEW/src/hypre/include -c -xSSE4.2 -O3 -D_LARGEFILE64_SOURCE -DMAXBLOCKS=2000 -DNXB=20 -DNYB=20 -DNZB=1 -DN_DIM=2", 400);
     strncpy(f_flags, "/nasa/sgi/mpt/2.15r20/bin/mpif90 -I/home3/svachhan/hdf5-1.8.20_NEW/hdf5/include -DH5_USE_16_API -I/home3/svachhan/hypre-2.11.2_NEW/src/hypre/include -c -r8 -i4 -real_size 64 -xSSE4.2 -align array32byte -O3 -DH5_USE_16_API -D_LARGEFILE64_SOURCE -D_FORTIFY_SOURCE=2 -DMAXBLOCKS=2000 -DNXB=20 -DNYB=20 -DNZB=1 -DN_DIM=2", 400);


}

