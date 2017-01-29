!!****h* source/physics/sourceTerms/Burn/Burn_interface
!!
!! This is the header file for the Burn module
!! that defines its public interfaces.
!!***

Module Burn_interface

#include "constants.h"
#include "Flash.h"

  interface Burn_computeDt
     subroutine Burn_computeDt (block_no,  &
          blkLimits,blkLimitsGC,  &
          solnData,   &
          dt_burn, dt_minloc )

       integer, intent(IN) :: block_no
       integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
       real, pointer, dimension(:,:,:,:) :: solnData
       real,INTENT(INOUT)    :: dt_burn
       integer,INTENT(INOUT) :: dt_minloc(5)
     end subroutine Burn_computeDt
  end interface

  interface Burn
     subroutine Burn(blockCount,blockList,dt)
       integer, intent(IN) :: blockCount
       integer,dimension(blockCount), intent(IN) :: blockList
       real,intent(IN) :: dt
     end subroutine Burn
  end interface

  interface Burn_finalize
     subroutine Burn_finalize()
     end subroutine Burn_finalize
  end interface

  interface Burn_init
     subroutine Burn_init()
       
     end subroutine Burn_init
  end interface

  interface Burn_nseAtDens
     subroutine Burn_nseAtDens(qbar_nse,sumyi_nse,approxtemp,edot,Yedot, Ye, dens, emq)
       implicit none
       real, intent(IN) :: Ye, dens, emq
       real, intent(OUT) :: qbar_nse,sumyi_nse,approxtemp,edot,Yedot
     end subroutine Burn_nseAtDens
  end interface


  interface Burn_nseAtPres
     subroutine Burn_nseAtPres(qbar_nse,sumyi_nse,approxtemp,edot,Yedot, Ye, pres, hmq)
       implicit none
       real, intent(IN)    :: Ye, pres, hmq
       real, intent(OUT)   :: qbar_nse,sumyi_nse,approxtemp,edot,Yedot
     end subroutine Burn_nseAtPres
  end interface


end Module Burn_interface


