! subroutine compute_errors
! Computes velocity and pressure error for Taylor Green Vortex problem.
!
! Marcos Vanella, July 2010
! --------------------------------------------------------------------------

  subroutine compute_errors(mype,blockCount,blockList)

  ! Modules Use:
  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkBC, &
    Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkIndexLimits, &
    Grid_getBlkBoundBox, Grid_getBlkCenterCoords

  use gr_interface, ONLY : gr_findMean

  use Driver_data, ONLY : dr_simTime

  use IncompNS_data, ONLY : ins_invRe
      
  implicit none
#include "constants.h"
#include "Flash.h"
#include "IncompNS.h"
  include 'mpif.h'

  integer :: mype
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(blockCount) :: blockList

  integer, parameter ::  ng = NGUARD
  integer, parameter ::  nxi= NGUARD + NXB
  integer, parameter ::  nyj= NGUARD + NYB
  integer, parameter ::  nzk= NGUARD + NZB
  integer, parameter ::  nxc= NGUARD + NXB + 1
  integer, parameter ::  nyc= NGUARD + NYB + 1
  integer, parameter ::  nzc= NGUARD + NZB + 1

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData
  
  real del(3),dx,dy,dz

  real meanPres

  integer lb,blockID,ierr

  real erru_inf,errv_inf,errp_inf,erru_L2,errv_L2,errp_L2
  integer Nux,Nvy,Npp
  real    :: mvisc

  real xedge,xcell,yedge,ycell

  integer i,j

  real coord(MDIM),bsize(MDIM)
  real, dimension(2,MDIM) :: boundBox

  ! -------------------------------------------------------------------

  mvisc = ins_invRe

  erru_inf =0.
  erru_L2  =0.
  Nux = 0
  errv_inf =0.
  errv_L2  =0.
  Nvy = 0
  errp_inf =0.
  errp_L2  =0.
  Npp = 0

  call gr_findMean(PRES_VAR,2,.false.,meanPres)


  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     ! Get Coord and Bsize for the block:
     ! Bounding box:
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)

     call Grid_getBlkCenterCoords(blockId,coord)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,solnData,CENTER)


     do j = NGUARD,NYB+2*NGUARD
        yedge = coord(JAXIS) - bsize(JAXIS)/2.0 +    &
                real(j - NGUARD - 1)*del(JAXIS)
        ycell = yedge + del(JAXIS)/2.0

        do i = NGUARD,NXB+2*NGUARD
           xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + &
                   real(i - NGUARD - 1)*del(IAXIS)
           xcell = xedge + del(IAXIS)/2.0

           facexData(RHDS_FACE_VAR,i,j,1) = -EXP(-2.0*mvisc*dr_simTime)* &
                                             COS(xedge)*SIN(ycell)

           faceyData(RHDS_FACE_VAR,i,j,1) =  EXP(-2.0*mvisc*dr_simTime)* &
                                             SIN(xcell)*COS(yedge)

           solnData(DELP_VAR,i,j,1)  = -0.25*EXP(-4.0*mvisc*dr_simTime)* &
                        ( COS(2.0*xcell) + COS(2.0*ycell) ) + meanPres 

        enddo
     enddo
     


     erru_L2  = erru_L2 + &
                sqrt(sum((facexData(VELC_FACE_VAR,ng+1:nxc,ng+1:nyj,1)-  &
                          facexData(RHDS_FACE_VAR,ng+1:nxc,ng+1:nyj,1))**2.))
     Nux = Nux + (NXB+1)*NYB


     erru_inf = max(erru_inf, &
                maxval(abs(facexData(VELC_FACE_VAR,ng+1:nxc,ng+1:nyj,1)- &
                           facexData(RHDS_FACE_VAR,ng+1:nxc,ng+1:nyj,1))))


     errv_L2  = errv_L2 + &
                sqrt(sum((faceyData(VELC_FACE_VAR,ng+1:nxi,ng+1:nyc,1)-  &
                          faceyData(RHDS_FACE_VAR,ng+1:nxi,ng+1:nyc,1))**2.))
     Nvy = Nvy + NXB*(NYB+1)


     errv_inf = max(errv_inf, &
                maxval(abs(faceyData(VELC_FACE_VAR,ng+1:nxi,ng+1:nyc,1)- &
                           faceyData(RHDS_FACE_VAR,ng+1:nxi,ng+1:nyc,1))))


     errp_L2  = errp_L2 + &
                sqrt(sum((solnData(PRES_VAR,ng+1:nxi,ng+1:nyj,1)- &
                          solnData(DELP_VAR,ng+1:nxi,ng+1:nyj,1))**2.))
     Npp = Npp + NXB*NYB


     errp_inf = max(errp_inf, &
                maxval(abs(solnData(PRES_VAR,ng+1:nxi,ng+1:nyj,1)- &
                           solnData(DELP_VAR,ng+1:nxi,ng+1:nyj,1))))



     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
           
  enddo

  if (mype .eq. 0) write(*,*) 'Nux,Nvy,Npp=',Nux,Nvy,Npp

  erru_L2  = erru_L2/sqrt(real(Nux))
  errv_L2  = errv_L2/sqrt(real(Nvy))
  errp_L2  = errp_L2/sqrt(real(Npp))


  if (mype .eq. 0) then

  !write(*,*) 'U Velocities:'
  write(*,*) mype,'einf U=',erru_inf,'eL2 U=',erru_L2
  !write(*,*) 'V Velocities:'
  write(*,*) mype,'einf V=',errv_inf,'eL2 U=',errv_L2
  !write(*,*) 'P Presures:'
  write(*,*) mype,'einf P=',errp_inf,'eL2 P=',errp_L2
  
  endif


  return
  end subroutine


