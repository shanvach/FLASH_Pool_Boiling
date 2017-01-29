subroutine ins_velomg2center(blockList,blockCount)


  use Grid_interface, only : Grid_getDeltas,         &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr


  implicit none
#include "constants.h"
#include "Flash.h"

  integer,intent(in) :: blockCount
  integer,intent(in) :: blockList(MAXBLOCKS)

  !Stub Function, to write cell centered velocities and vorticity velx,vely,velz,omgx,omgy,omgz
  ! add REQUIRES physics/IncompNS/IncompNSMain/extras
  ! These will be writen to the flash plotfiles if required in
  ! flash.par . 

  return

end subroutine ins_velomg2center
