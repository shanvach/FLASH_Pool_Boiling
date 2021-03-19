! Subroutine outtotecplot
!
! Subroutine to write out to Tecplot data in binary form.
!
! ---------------------------------------------------------------------------

#include "constants.h"
#include "Flash.h"

      subroutine outtotecplot_uv(mype,time,dt,istep,count, &
                 timer,blockList,blockCount,firstfileflag)

      use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkPtr, &
        Grid_releaseBlkPtr, Grid_getBlkIndexLimits, &
        Grid_getBlkBoundBox, Grid_getBlkCenterCoords


#ifdef FLASH_GRID_UG
#else
      use physicaldata, ONLY : interp_mask_unk,interp_mask_unk_res
#endif


  use IncompNS_data, ONLY : ins_invRe


  implicit none

  include "Flash_mpi.h"
      
  integer, intent(in) :: mype,istep,count,firstfileflag
  integer, intent(in) :: blockCount
  integer, intent(in) :: blockList(MAXBLOCKS)
  real, intent(in) :: time,dt,timer

      
  ! Local Variables:
  integer :: numblocks,var,i,j,k,lb,nxc,nyc,nzc
  character(25) :: filename
  character(6) :: index_lb,index_mype

  real xedge(NXB+1),xcell(NXB+1)
  real yedge(NYB+1),ycell(NYB+1)
  real intsx(NXB+1), intsy(NYB+1)


  real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData

  real facevarxx(NXB+2*NGUARD+1,NYB+2*NGUARD),   &
       facevarxxan(NXB+2*NGUARD+1,NYB+2*NGUARD), &
       difffacevarx(NXB+2*NGUARD+1,NYB+2*NGUARD)

  real facevaryy(NXB+2*NGUARD,NYB+2*NGUARD+1),   &
       facevaryyan(NXB+2*NGUARD,NYB+2*NGUARD+1), & 
       difffacevary(NXB+2*NGUARD,NYB+2*NGUARD+1)

  real phivar(NXB+2*NGUARD,NYB+2*NGUARD)

  real presvar(NXB+2*NGUARD,NYB+2*NGUARD),       &
       presvar_an(NXB+2*NGUARD,NYB+2*NGUARD),    &
       diffpresvar(NXB+2*NGUARD,NYB+2*NGUARD)

  real xi,yj

  real*4 arraylb(NXB+1,NYB+2,1),        &
         arraylb2(NXB+2,NYB+1,1),       &
         arraylb3(NXB,NYB,1)

  integer blockID
  real del(MDIM),dx,dy
  real, dimension(MDIM)  :: coord,bsize
  real ::  boundBox(2,MDIM)

  integer*4 TecIni,TecDat,TecZne,TecNod,TecFil,TecEnd
  integer*4 VIsdouble
  integer*4 Debug,ijk,ijk2,ijk3,Npts,NElm
  character*1 NULLCHR
  
  real :: pi, tn, mvisc

    end subroutine outtotecplot_uv 
