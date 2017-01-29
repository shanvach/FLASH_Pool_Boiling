subroutine setptindices(opind,npind,ovind,nvind,intvind,ofind,nfind)

  implicit none
#include "Flash.h"
#include "constants.h"
  
  integer,dimension(3),INTENT(out)::opind,npind,ovind,nvind,intvind,ofind,nfind
  
 ! Set the indices for variables at n and n+1
  ! Old positions and velocities indices
  opind(IAXIS)= POSX_PART_PROP
  opind(JAXIS)= POSY_PART_PROP
  opind(KAXIS)= POSZ_PART_PROP
  ovind(IAXIS)= VELX_PART_PROP
  ovind(JAXIS)= VELY_PART_PROP
  ovind(KAXIS)= VELZ_PART_PROP
  ofind(IAXIS)= FNX_PART_PROP
  ofind(JAXIS)= FNY_PART_PROP
  ofind(KAXIS)= FNZ_PART_PROP  
  ! N positions and velocities indices
  npind(IAXIS)= PNPX_PART_PROP
  npind(JAXIS)= PNPY_PART_PROP
  npind(KAXIS)= PNPZ_PART_PROP
  nvind(IAXIS)= VNPX_PART_PROP
  nvind(JAXIS)= VNPY_PART_PROP
  nvind(KAXIS)= VNPZ_PART_PROP
  nfind(IAXIS)= FNPX_PART_PROP
  nfind(JAXIS)= FNPY_PART_PROP
  nfind(KAXIS)= FNPZ_PART_PROP
  intvind(IAXIS)=VIX_PART_PROP
  intvind(JAXIS)=VIY_PART_PROP
  intvind(KAXIS)=VIZ_PART_PROP
  

end subroutine setptindices
