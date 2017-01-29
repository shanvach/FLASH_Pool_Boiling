!!****if* source/Simulation/SimulationMain/INavierStokes/3D/ChannelLam/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!!  Reference:
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!  myPE   -           my processor number
!!
!! 
!!
!!***

subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY : sim_xMin, sim_xMax, &
                              sim_yMin, sim_yMax, &
                              sim_zMin, sim_zMax, &
                              sim_gCell

  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkBoundBox,    &
                             Grid_getBlkCenterCoords

  use Grid_data, only : gr_meshMe !,gr_axisComm,gr_axisMe,gr_axisNumProcs,gr_imax,gr_imin

  use Driver_interface, only : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  integer, intent(in) :: blockID
  !!$ ---------------------------------
 
  integer :: i, j, k
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) ::  blIndSize,blIndSizeGC

  real, dimension(MDIM)  :: coord,bsize,del
  real ::  boundBox(2,MDIM)
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData


  integer :: nxvar,nyvar,nzvar,nxvar2,nyvar2,nzvar2,auxint,var
  integer :: npt
  real, allocatable, dimension(:)   :: xvar,yvar,zvar,aux
  real, allocatable, dimension(:)   :: xcen,ycen,zcen
  real, allocatable, dimension(:,:,:) :: Var0u,Var0v,Var0w
  real, allocatable, dimension(:,:,:) :: Var0p,uab1,vab1,wab1

  real :: texternal,qo,qold,dpdx

  real :: dx,dy,dz
  integer :: ix1,ix2,jy1,jy2,kz1,kz2

  character(20) :: filename

  real :: xedge(GRID_IHI_GC+1),xcell(GRID_IHI_GC)
  real :: yedge(GRID_JHI_GC+1),ycell(GRID_JHI_GC)
  real :: zedge(GRID_KHI_GC+1),zcell(GRID_KHI_GC)

  real :: zp, fb

  ! ISOTURB_UG_5PI_512
  real, parameter :: ao = 0.1
  real, parameter :: distb = 3.*PI/2.

  ! ISOTURB_UG_4PI_TEST
  !real, parameter :: ao = 0.1
  !real, parameter :: distb = 6.*3.*PI/2. !OJO MADE IT NINE PI TO COVER ALL BOX

  ! ISOTURB_UG_5PI_4OVER25_512
  !real, parameter :: ao = 0.1*4./25. ! Rescale constant through (2pi)^2/(5pi)^2
  !real, parameter :: distb = 3.*PI/2.

  ! ISOTURB_UG_5PI_2OVER5_512
  !real, parameter :: ao = 0.1*2./5. ! Rescale constant through (2pi)/(5pi)
  !real, parameter :: distb = 3.*PI/2.


  real, parameter :: distd = PI/2.

  real, parameter :: zforcecen = 0.

  real, parameter :: zFreeSurface = 1. * 2.5*PI

  !----------------------------------------------------------------------

  ! Get Coord and Bsize for the block:
  ! Bounding box:
  call Grid_getBlkBoundBox(blockId,boundBox)
  bsize(:) = boundBox(2,:) - boundBox(1,:)

  call Grid_getBlkCenterCoords(blockId,coord)

  ! Get blocks dx, dy ,dz:
  call Grid_getDeltas(blockID,del)

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)

  ! Compute Grid line locations:
  ! Z:
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     zedge(k) = coord(KAXIS) - bsize(KAXIS)/2.0 + real(k - NGUARD - 1)*del(KAXIS)
     zcell(k) = zedge(k) + 0.5*del(KAXIS)
  enddo
  zedge(blkLimitsGC(HIGH,KAXIS)+1)=coord(KAXIS)-bsize(KAXIS)/2.0 + &
                                   real(blkLimitsGC(HIGH,KAXIS)-NGUARD)*del(KAXIS)
  ! Y:
  do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
     yedge(j) = coord(JAXIS) - bsize(JAXIS)/2.0 + real(j - NGUARD - 1)*del(JAXIS)
     ycell(j) = yedge(j) + 0.5*del(JAXIS)
  enddo
  yedge(blkLimitsGC(HIGH,JAXIS)+1)=coord(JAXIS)-bsize(JAXIS)/2.0 + &
                                   real(blkLimitsGC(HIGH,JAXIS)-NGUARD)*del(JAXIS)

  ! X:
  do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
     xedge(i) = coord(IAXIS) - bsize(IAXIS)/2.0 + real(i - NGUARD - 1)*del(IAXIS)
     xcell(i) = xedge(i) + 0.5*del(IAXIS)
  enddo
  xedge(blkLimitsGC(HIGH,IAXIS)+1)=coord(IAXIS)-bsize(IAXIS)/2.0+ &
                                   real(blkLimitsGC(HIGH,IAXIS)-NGUARD)*del(IAXIS)


  ! Point to Blocks centered variables:
  call Grid_getBlkPtr(blockID,solnData,CENTER)

  ! Point to Blocks face variables: 
  call Grid_getBlkPtr(blockID,facexData,FACEX)
  call Grid_getBlkPtr(blockID,faceyData,FACEY)
  call Grid_getBlkPtr(blockID,facezData,FACEZ)

  ! Initialize Distance Function:
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

          solnData(DFUN_VAR,i,j,k) = zcell(k) - zFreeSurface
 
        enddo
     enddo
  enddo


  ! Allocate arrays to read variable:
  npt = 64
  nxvar = npt+2
  nyvar = npt+2
  nzvar = npt+2

  allocate(xvar(nxvar),yvar(nyvar),zvar(nzvar),aux(nzvar))
  allocate(xcen(nxvar),ycen(nyvar),zcen(nzvar))
  allocate(Var0u(nxvar,nyvar,nzvar), &
           Var0v(nxvar,nyvar,nzvar), &
           Var0w(nxvar,nyvar,nzvar), &
           Var0p(nxvar,nyvar,nzvar), &
           uab1(nxvar,nyvar,nzvar),  &
           vab1(nxvar,nyvar,nzvar),  &
           wab1(nxvar,nyvar,nzvar))

  ! read grid:
  OPEN(UNIT=2,FILE='./ICdata/64x64x64.grd',STATUS='OLD' )
  read(2,111) ix1,ix2
  read(2,111) jy1,jy2
  read(2,111) kz1,kz2
  read(2,112) dx,dy,dz

  read(2,*) (xvar(i),i=ix1,ix2+1)
  read(2,*) (aux(i) ,i=ix1,ix2+1)
  read(2,*) (aux(i) ,i=ix1,ix2+1)
  read(2,*) (aux(i) ,i=ix1,ix2+1)
  read(2,*) (aux(i) ,i=ix1,ix2+1)
  read(2,*) (aux(i) ,i=ix1,ix2+1)

  read(2,*) (yvar(j),j=jy1,jy2+1)
  read(2,*) (aux(j) ,j=jy1,jy2+1)
  read(2,*) (aux(j) ,j=jy1,jy2+1)
  read(2,*) (aux(j) ,j=jy1,jy2+1)
  read(2,*) (aux(j) ,j=jy1,jy2+1)
  read(2,*) (aux(j) ,j=jy1,jy2+1)

  read(2,*) (zvar(k) ,k=kz1,kz2+1)
  read(2,*) (aux(k),k=kz1,kz2+1)
  read(2,*) (aux(k) ,k=kz1,kz2+1)
  read(2,*) (aux(k) ,k=kz1,kz2+1)
  read(2,*) (aux(k) ,k=kz1,kz2+1)
  read(2,*) (aux(k) ,k=kz1,kz2+1)

  close(2)

  xvar(npt+1) = 2.*xvar(npt) - xvar(npt-1)
  yvar(npt+1) = 2.*yvar(npt) - yvar(npt-1)
  zvar(npt+1) = 2.*zvar(npt) - zvar(npt-1)

  ! Compute cell centered coordinates:
  do i = 2,ix2
     xcen(i) = 0.5*(xvar(i)+xvar(i-1))
  enddo
  xcen(1) = xcen(2) - (xcen(3)-xcen(2))
  xcen(nxvar) = xcen(nxvar-1) + (xcen(nxvar-1)-xcen(nxvar-2))

  do j = 2,jy2
     ycen(j) = 0.5*(yvar(j)+yvar(j-1))
  enddo
  ycen(1) = ycen(2) - (ycen(3)-ycen(2))
  ycen(nyvar) = ycen(nyvar-1) + (ycen(nyvar-1)-ycen(nyvar-2))

  do k = 2,kz2
     zcen(k) = 0.5*(zvar(k)+zvar(k-1))
  enddo
  zcen(1) = zcen(2) - (zcen(3)-zcen(2))
  zcen(nzvar) = zcen(nzvar-1) + (zcen(nzvar-1)-zcen(nzvar-2))

  ! read restart file
  open(19,file='./ICdata/64x64x64.res', status='unknown',form='unformatted')
  read(19)  texternal
  read(19)  qo,qold,dpdx
  read(19)  Var0u,Var0v,Var0w,Var0p
  read(19)  uab1,vab1,wab1
  close(19)


  ! Substract pi to match new box position:
  xvar = xvar - 4.*atan(1.)
  yvar = yvar - 4.*atan(1.)
  zvar = zvar - 4.*atan(1.)

  xcen = xcen - 4.*atan(1.)
  ycen = ycen - 4.*atan(1.)
  zcen = zcen - 4.*atan(1.)

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,FACEX)
  ! Interpolate the variable to each point in the AMR grid
  call interpvars_face(NGUARD,NXB,NYB,NZB, &
                     coord,bsize,&
                     nxvar,nyvar,nzvar,xvar,ycen,zcen,&
                     Var0u,blkLimitsGC(HIGH,IAXIS),&
                     blkLimitsGC(HIGH,JAXIS),&
                     blkLimitsGC(HIGH,KAXIS),&
                     facexData(VELC_FACE_VAR,:,:,:),1)


  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,FACEY)
  ! Interpolate the variable to each point in the AMR grid
  call interpvars_face(NGUARD,NXB,NYB,NZB, &
                     coord,bsize,&
                     nxvar,nyvar,nzvar,xcen,yvar,zcen,&
                     Var0v,blkLimitsGC(HIGH,IAXIS),&
                     blkLimitsGC(HIGH,JAXIS),&
                     blkLimitsGC(HIGH,KAXIS),&
                     faceyData(VELC_FACE_VAR,:,:,:),2)


  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,FACEZ)
  ! Interpolate the variable to each point in the AMR grid
  call interpvars_face(NGUARD,NXB,NYB,NZB,&
                      coord,bsize,&
                      nxvar,nyvar,nzvar,xcen,ycen,zvar,&
                      Var0w,blkLimitsGC(HIGH,IAXIS),&
                      blkLimitsGC(HIGH,JAXIS),&
                      blkLimitsGC(HIGH,KAXIS),&
                      facezData(VELC_FACE_VAR,:,:,:),3)


  ! deallocation of auxiliary arrays exit
  deallocate(xvar,yvar,zvar,aux)
  deallocate(xcen,ycen,zcen)
  deallocate(Var0u,Var0v,Var0w,Var0p)
  deallocate(uab1,vab1,wab1)

!  write(*,*) 'BlockID=',blockID
!  write(*,*) 'Center coordinates=',coord
!  write(*,*) 'Size =',bsize
!  write(*,*) 'Uveloc, 2,2,2=', &
!              facexData(VELC_FACE_VAR,NGUARD+1+2,NGUARD+2,NGUARD+2)
!  write(*,*) 'Vveloc, 2,2,2=', &
!              faceyData(VELC_FACE_VAR,NGUARD+2,NGUARD+1+2,NGUARD+2)
!  write(*,*) 'Vveloc, 2,2,2=', &
!              facezData(VELC_FACE_VAR,NGUARD+2,NGUARD+2,NGUARD+1+2)


  solnData(PRES_VAR,:,:,:) = 0.0
  solnData(DELP_VAR,:,:,:) = 0.0
  solnData(DUST_VAR,:,:,:) = 0.0
  solnData(TVIS_VAR,:,:,:) = 0.0

  facexData(RHDS_FACE_VAR,:,:,:) = 0.0
  faceyData(RHDS_FACE_VAR,:,:,:) = 0.0
  facezData(RHDS_FACE_VAR,:,:,:) = 0.0

  ! Set Force distribution Fb(x,y,z):
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)
  ! X velocity:
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1

          zp = abs(zcell(k)-zforcecen) 

          if (zp .le. distb) then
             fb = 1.
          elseif (zp .le. (distb+distd)) then
             fb = 0.5*(1.-cos(PI/distd*(zp-distb-distd)))
          else
             fb = 0.
          endif

          facexData(FBAO_FACE_VAR,i,j,k) = ao*fb
        enddo
     enddo
  enddo

  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

          zp = abs(zcell(k)-zforcecen)

          if (zp .le. distb) then
             fb = 1.
          elseif (zp .le. (distb+distd)) then
             fb = 0.5*(1.-cos(PI/distd*(zp-distb-distd)))
          else
             fb = 0.
          endif
       
          faceyData(FBAO_FACE_VAR,i,j,k) = ao*fb         

        enddo
     enddo
  enddo

  ! Z velocity:
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)+1
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

          zp = abs(zedge(k)-zforcecen)

          if (zp .le. distb) then
             fb = 1.
          elseif (zp .le. (distb+distd)) then
             fb = 0.5*(1.-cos(PI/distd*(zp-distb-distd)))
          else
             fb = 0.
          endif

          facezData(FBAO_FACE_VAR,i,j,k) = ao*fb

        enddo
     enddo
  enddo


  ! Divergence:
  solnData(DUST_VAR,:,:,:) = 0.0
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
          solnData(DUST_VAR,i,j,k) =(facexData(VELC_FACE_VAR,i+1,j,k)-facexData(VELC_FACE_VAR,i,j,k))/del(IAXIS)+&
                                    (faceyData(VELC_FACE_VAR,i,j+1,k)-faceyData(VELC_FACE_VAR,i,j,k))/del(JAXIS)+&
                                    (facezData(VELC_FACE_VAR,i,j,k+1)-facezData(VELC_FACE_VAR,i,j,k))/del(KAXIS)

        enddo
     enddo
  enddo
  if (gr_meshMe .eq. MASTER_PE)      write(*,*) 'Mype=',gr_meshMe,', min,maxdiv=', &
  minval(solnData(DUST_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)),&
  maxval(solnData(DUST_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI))


  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
  call Grid_releaseBlkPtr(blockID,facezData,FACEZ)


  return

111    format (i4,3x,i4)
112    format (3(3x,e12.4))


end subroutine Simulation_initBlock
