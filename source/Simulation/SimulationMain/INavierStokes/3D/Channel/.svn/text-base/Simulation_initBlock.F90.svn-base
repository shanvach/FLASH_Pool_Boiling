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

!!#define LAMINAR_FLOW
#define DOUBLE_GRID

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

  use Grid_data, only : gr_meshMe,gr_axisComm,gr_axisMe,gr_axisNumProcs,gr_imax,gr_imin

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

  real :: xedge(GRID_IHI_GC+1),xcell(GRID_IHI_GC)
  real :: yedge(GRID_JHI_GC+1),ycell(GRID_JHI_GC)
  real :: zedge(GRID_KHI_GC+1),zcell(GRID_KHI_GC)

  real, parameter :: eps = 1.e-14

  real :: du,dv,dw
 
  real, parameter :: Uoo = 23.3333 ! Laminar umax for ubulk 15.55, channel height goes
                                   ! from -1 to 1.

  real :: wx,wz,ampi,kwi,factor
  integer :: ik

  integer, parameter :: kwlen = 2

  real :: kw(kwlen),amp(kwlen)

  real, allocatable, dimension(:,:,:) :: ucoarse,vcoarse,wcoarse

  integer :: kstart,factorz,nprocy_coarse,nprocz_coarse,iprocy_coarse,iprocz_coarse,proc_coarse

  integer :: NXBC , NYBC , NZBC , grd_ihi , grd_jhi , grd_khi , grd_ihi_gc , grd_jhi_gc , grd_khi_gc     

  character(20) :: filename

  !----------------------------------------------------------------------



  kw(1) = 2.; kw(2) = 4.;
  amp(1)= 4.; amp(2)= 2.;

  ! Get Coord and Bsize for the block:
  ! Bounding box:
  call Grid_getBlkBoundBox(blockId,boundBox)
  bsize(:) = boundBox(2,:) - boundBox(1,:)

  call Grid_getBlkCenterCoords(blockId,coord)

  ! Get blocks dx, dy ,dz:
  call Grid_getDeltas(blockID,del)

  factor = del(IAXIS)/del(KAXIS)
  if (gr_meshMe .eq. MASTER_PE) write(*,*) 'DELZ=',del(KAXIS),', DELX=',del(IAXIS),', DELZ/DELX=',del(KAXIS)/del(IAXIS)


  ! Point to Blocks centered variables:
  call Grid_getBlkPtr(blockID,solnData,CENTER)

  ! Point to Blocks face variables: 
  call Grid_getBlkPtr(blockID,facexData,FACEX)
  call Grid_getBlkPtr(blockID,faceyData,FACEY)
  call Grid_getBlkPtr(blockID,facezData,FACEZ)

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

  facexData(VELC_FACE_VAR,:,:,:) = 0.0
  faceyData(VELC_FACE_VAR,:,:,:) = 0.0
  facezData(VELC_FACE_VAR,:,:,:) = 0.0

#ifdef LAMINAR_FLOW

  !CALL RANDOM_SEED

  !! Velocities: 
  !! X velocity:
  !do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
  !   do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
  !      do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1
  !          call random_number(du)
  !          facexData(VELC_FACE_VAR,i,j,k) = eps*du
  !      enddo
  !   enddo
  !enddo

  !! X velocity:
  !do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
  !   do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
  !      do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
  !          call random_number(dv)
  !          faceyData(VELC_FACE_VAR,i,j,k) = eps*dv
  !      enddo
  !   enddo
  !enddo

  ! Z velocity:
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)+1
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
            !call random_number(dw)            
            facezData(VELC_FACE_VAR,i,j,k) = Uoo*(1-ycell(j)**2.) !+ eps*dw
        enddo
     enddo
  enddo

  do ik=1,kwlen

    kwi =kw(ik)
    ampi=amp(ik)

    wx = 2*PI*kwi/(gr_imax-gr_imin)
    wz = factor*wx

    ! X velocity:
    do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
       do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
          do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1
            du = -ampi*cos(wz*zcell(k))*sin(wx*xedge(i))
            facexData(VELC_FACE_VAR,i,j,k) = facexData(VELC_FACE_VAR,i,j,k) + du
          enddo
       enddo
    enddo

    ! Z velocity:
    do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)+1
       do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
          do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
            dw = ampi/factor*sin(wz*zedge(k))*cos(wx*xcell(i))
            facezData(VELC_FACE_VAR,i,j,k) = facezData(VELC_FACE_VAR,i,j,k) + dw
          enddo
       enddo
    enddo

  enddo

  ! Divergence:
  solnData(DUST_VAR,:,:,:) = 0.0
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
          solnData(DUST_VAR,i,j,k) = (facexData(VELC_FACE_VAR,i+1,j,k)-facexData(VELC_FACE_VAR,i,j,k))/del(IAXIS) + &
                                     (faceyData(VELC_FACE_VAR,i,j+1,k)-faceyData(VELC_FACE_VAR,i,j,k))/del(JAXIS) + &
                                     (facezData(VELC_FACE_VAR,i,j,k+1)-facezData(VELC_FACE_VAR,i,j,k))/del(KAXIS)
        enddo
     enddo
  enddo
  write(*,*) 'Mype=',gr_meshMe,', min,max div=',minval(solnData(DUST_VAR,:,:,:)),maxval(solnData(DUST_VAR,:,:,:))

#else

#ifdef DOUBLE_GRID

  NXBC = NXB/2
  NYBC = NYB/2
  NZBC = NZB/2

  grd_ihi    = NXBC + NGUARD
  grd_jhi    = NYBC + NGUARD
  grd_khi    = NZBC + NGUARD

  grd_ihi_gc = NXBC + 2*NGUARD
  grd_jhi_gc = NYBC + 2*NGUARD
  grd_khi_gc = NZBC + 2*NGUARD

  allocate(ucoarse(GRID_ILO_GC:grd_ihi_gc+1,GRID_JLO_GC:grd_jhi_gc,GRID_KLO_GC:grd_khi_gc))
  allocate(vcoarse(GRID_ILO_GC:grd_ihi_gc,GRID_JLO_GC:grd_jhi_gc+1,GRID_KLO_GC:grd_khi_gc))
  allocate(wcoarse(GRID_ILO_GC:grd_ihi_gc,GRID_JLO_GC:grd_jhi_gc,GRID_KLO_GC:grd_khi_gc+1))

  ! Same number of processors from coarse read and finer grid
  nprocy_coarse = gr_axisNumProcs(JAXIS)
  nprocz_coarse = gr_axisNumProcs(KAXIS)

  iprocy_coarse = gr_axisMe(JAXIS)
  iprocz_coarse = gr_axisMe(KAXIS)

  proc_coarse = iprocz_coarse*nprocy_coarse + iprocy_coarse 

  if (proc_coarse .ne. gr_meshMe) then
    write(*,*) 'gr_meshMe=',gr_meshMe,', proc_coarse=',proc_coarse
    call Driver_abortFlash("Proc_coarse ne gr_meshMe.") 
  endif

  ! file:
  write(filename, '("IOData/Restrt.", i5.5)') proc_coarse
  open(33,file=trim(filename),form='unformatted')

  ! U velocity:
  do k=GRID_KLO,grd_khi
    do j=GRID_JLO,grd_jhi
      do i=GRID_ILO,grd_ihi+1

        read(33) ucoarse(i,j,k)

      enddo
    enddo
  enddo

  ! V velocity:
  do k=GRID_KLO,grd_khi
    do j=GRID_JLO,grd_jhi+1
      do i=GRID_ILO,grd_ihi

        read(33) vcoarse(i,j,k)

      enddo
    enddo
  enddo

  ! W velocity:
  do k=GRID_KLO,grd_khi+1
    do j=GRID_JLO,grd_jhi
      do i=GRID_ILO,grd_ihi

        read(33) wcoarse(i,j,k)

      enddo
    enddo
  enddo

  close(33)

  if( gr_meshMe .eq. MASTER_PE) write(*,*) 'Read Binary Start Files.'

  ! Set faces:
  ! U velocity:
  do k=GRID_KLO,GRID_KHI
    do j=GRID_JLO,GRID_JHI
      do i=GRID_ILO,GRID_IHI+2

        facexData(VELC_FACE_VAR,i,j,k) = ucoarse(GRID_ILO+(i-GRID_ILO)/2, &
                                                 GRID_JLO+(j-GRID_JLO)/2, &
                                                 GRID_KLO+(k-GRID_KLO)/2)

      enddo
    enddo
  enddo

  ! V velocity:
  do k=GRID_KLO,GRID_KHI
    do j=GRID_JLO,GRID_JHI+2
      do i=GRID_ILO,GRID_IHI

        faceyData(VELC_FACE_VAR,i,j,k) = vcoarse(GRID_ILO+(i-GRID_ILO)/2, &
                                                 GRID_JLO+(j-GRID_JLO)/2, &
                                                 GRID_KLO+(k-GRID_KLO)/2)

      enddo
    enddo
  enddo
  
  ! W velocity:
  do k=GRID_KLO,GRID_KHI+2
    do j=GRID_JLO,GRID_JHI
      do i=GRID_ILO,GRID_IHI

        facezData(VELC_FACE_VAR,i,j,k) = wcoarse(GRID_ILO+(i-GRID_ILO)/2, &
                                                 GRID_JLO+(j-GRID_JLO)/2, &
                                                 GRID_KLO+(k-GRID_KLO)/2)

      enddo
    enddo
  enddo

  ! ReWrite internal values:
  ! U velocities:
  do k=GRID_KLO,GRID_KHI
    do j=GRID_JLO,GRID_JHI
      do i=GRID_ILO+1,GRID_IHI,2

        facexData(VELC_FACE_VAR,i,j,k) = 0.5*(facexData(VELC_FACE_VAR,i-1,j,k)+facexData(VELC_FACE_VAR,i+1,j,k))

      enddo
    enddo
  enddo
  
  ! V velocities:
  do k=GRID_KLO,GRID_KHI
    do j=GRID_JLO+1,GRID_JHI,2
      do i=GRID_ILO,GRID_IHI

        faceyData(VELC_FACE_VAR,i,j,k) = 0.5*(faceyData(VELC_FACE_VAR,i,j-1,k)+faceyData(VELC_FACE_VAR,i,j+1,k))

      enddo
    enddo
  enddo
  
  ! W velocities:
  do k=GRID_KLO+1,GRID_KHI,2
    do j=GRID_JLO,GRID_JHI
      do i=GRID_ILO,GRID_IHI

        facezData(VELC_FACE_VAR,i,j,k) = 0.5*(facezData(VELC_FACE_VAR,i,j,k-1)+facezData(VELC_FACE_VAR,i,j,k+1))

      enddo
    enddo  
  enddo

#else
  ! Read initial field from binary file and interpolate:
  factorz       = 2 
  nprocy_coarse = gr_axisNumProcs(JAXIS)
  nprocz_coarse = factorz*gr_axisNumProcs(KAXIS)

  iprocy_coarse = gr_axisMe(JAXIS)
  iprocz_coarse = gr_axisMe(KAXIS)/factorz

  proc_coarse   = iprocz_coarse*nprocy_coarse + iprocy_coarse 

  allocate(ucoarse(GRID_ILO_GC:GRID_IHI_GC+1,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC))
  allocate(vcoarse(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC+1,GRID_KLO_GC:GRID_KHI_GC))
  allocate(wcoarse(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC+1))


!  write(*,*) &
!'Mype=',gr_meshMe,gr_axisMe(JAXIS),gr_axisMe(KAXIS),iprocy_coarse,iprocz_coarse,proc_coarse

  ! file:
  write(filename, '("IOData/Restrt.", i5.5)') proc_coarse
  open(33,file=trim(filename),form='unformatted')

  ! U velocity:
  do k=GRID_KLO,GRID_KHI
    do j=GRID_JLO,GRID_JHI
      do i=GRID_ILO,GRID_IHI+1

        read(33) ucoarse(i,j,k)

      enddo
    enddo
  enddo

  ! V velocity:
  do k=GRID_KLO,GRID_KHI
    do j=GRID_JLO,GRID_JHI+1
      do i=GRID_ILO,GRID_IHI

        read(33) vcoarse(i,j,k)

      enddo
    enddo
  enddo

  ! W velocity:
  do k=GRID_KLO,GRID_KHI+1
    do j=GRID_JLO,GRID_JHI
      do i=GRID_ILO,GRID_IHI

        read(33) wcoarse(i,j,k)

      enddo
    enddo
  enddo

  close(33)

  ! Assign to finer cells:
  kstart = GRID_KLO
  if (mod(gr_axisMe(KAXIS)+1,factorz) .eq. 0) kstart = GRID_KLO + NZB/factorz

!  write(*,*) 'Mype=',gr_meshMe,gr_axisMe(JAXIS),gr_axisMe(KAXIS),kstart

  if (mod(NZB,2) .ne. 0) then

  if (kstart .eq. GRID_KLO) then
  ! Assign U velocities
  do k=GRID_KLO,GRID_KHI
    do j=GRID_JLO,GRID_JHI
       do i=GRID_ILO,GRID_IHI+1
         facexData(VELC_FACE_VAR,i,j,k) = ucoarse(i,j,kstart+(k-NGUARD-1)/factorz)
       enddo
    enddo
  enddo
  ! Assign V velocities
  do k=GRID_KLO,GRID_KHI
    do j=GRID_JLO,GRID_JHI+1
       do i=GRID_ILO,GRID_IHI
         faceyData(VELC_FACE_VAR,i,j,k) = vcoarse(i,j,kstart+(k-NGUARD-1)/factorz)
       enddo
    enddo
  enddo
  ! Assign W velocities
  do k=GRID_KLO-1,GRID_KHI+2
    do j=GRID_JLO,GRID_JHI
       do i=GRID_ILO,GRID_IHI
         facezData(VELC_FACE_VAR,i,j,k) = wcoarse(i,j,kstart+(k-NGUARD-1)/factorz)
       enddo
    enddo
  enddo
  ! Fix fine grid mid vels:
  do k=GRID_KLO+1,GRID_KHI+1,factorz
    do j=GRID_JLO,GRID_JHI
       do i=GRID_ILO,GRID_IHI
         facezData(VELC_FACE_VAR,i,j,k) = (facezData(VELC_FACE_VAR,i,j,k-1)+facezData(VELC_FACE_VAR,i,j,k+1))/real(factorz)
       enddo
    enddo
  enddo

  else  ! kstart = GRID_KLO + NZB/factorz

  ! Assign U velocities
  do k=GRID_KLO,GRID_KHI
    do j=GRID_JLO,GRID_JHI
       do i=GRID_ILO,GRID_IHI+1
         facexData(VELC_FACE_VAR,i,j,k) = ucoarse(i,j,kstart+(k-NGUARD)/factorz)
       enddo
    enddo
  enddo
  ! Assign V velocities
  do k=GRID_KLO,GRID_KHI
    do j=GRID_JLO,GRID_JHI+1
       do i=GRID_ILO,GRID_IHI
         faceyData(VELC_FACE_VAR,i,j,k) =vcoarse(i,j,kstart+(k-NGUARD)/factorz)
       enddo
    enddo
  enddo
  ! Assign W velocities
  do k=GRID_KLO-1,GRID_KHI+2
    do j=GRID_JLO,GRID_JHI
       do i=GRID_ILO,GRID_IHI
         facezData(VELC_FACE_VAR,i,j,k) =wcoarse(i,j,kstart+(k-NGUARD)/factorz)
       enddo
    enddo
  enddo
  ! Fix fine grid mid vels:
  do k=GRID_KLO,GRID_KHI,factorz
    do j=GRID_JLO,GRID_JHI
       do i=GRID_ILO,GRID_IHI
         facezData(VELC_FACE_VAR,i,j,k) =(facezData(VELC_FACE_VAR,i,j,k-1)+facezData(VELC_FACE_VAR,i,j,k+1))/real(factorz)
       enddo
    enddo
  enddo


  endif

  else

  ! Assign U velocities
  do k=GRID_KLO,GRID_KHI
    do j=GRID_JLO,GRID_JHI
       do i=GRID_ILO,GRID_IHI+1
         facexData(VELC_FACE_VAR,i,j,k) = ucoarse(i,j,kstart+(k-NGUARD-1)/factorz)
       enddo
    enddo
  enddo
  ! Assign V velocities
  do k=GRID_KLO,GRID_KHI
    do j=GRID_JLO,GRID_JHI+1
       do i=GRID_ILO,GRID_IHI
         faceyData(VELC_FACE_VAR,i,j,k) = vcoarse(i,j,kstart+(k-NGUARD-1)/factorz)
       enddo
    enddo
  enddo
  ! Assign W velocities
  do k=GRID_KLO-1,GRID_KHI+1
    do j=GRID_JLO,GRID_JHI
       do i=GRID_ILO,GRID_IHI
         facezData(VELC_FACE_VAR,i,j,k) = wcoarse(i,j,kstart+(k-NGUARD-1)/factorz)
       enddo
    enddo
  enddo
  ! Fix fine grid mid vels:
  do k=GRID_KLO+1,GRID_KHI,factorz
    do j=GRID_JLO,GRID_JHI
       do i=GRID_ILO,GRID_IHI
         facezData(VELC_FACE_VAR,i,j,k) =(facezData(VELC_FACE_VAR,i,j,k-1)+facezData(VELC_FACE_VAR,i,j,k+1))/real(factorz)
       enddo
    enddo
  enddo

  endif

#endif

  deallocate(ucoarse,vcoarse,wcoarse)

  ! Divergence:
  solnData(DUST_VAR,:,:,:) = 0.0
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
          solnData(DUST_VAR,i,j,k) = (facexData(VELC_FACE_VAR,i+1,j,k)-facexData(VELC_FACE_VAR,i,j,k))/del(IAXIS) + &
                                     (faceyData(VELC_FACE_VAR,i,j+1,k)-faceyData(VELC_FACE_VAR,i,j,k))/del(JAXIS) + &
                                     (facezData(VELC_FACE_VAR,i,j,k+1)-facezData(VELC_FACE_VAR,i,j,k))/del(KAXIS)

        enddo
     enddo
  enddo
  write(*,*) 'Mype=',gr_meshMe,', min,maxdiv=', &
minval(solnData(DUST_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)), &
maxval(solnData(DUST_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI))

! if(gr_meshMe .eq. MASTER_PE) then
!  write(*,*) 'Mype=',gr_meshMe,', min,maxdiv=', &
!minval(solnData(DUST_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)),&
!maxval(solnData(DUST_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI))
!  do k=GRID_KLO,GRID_KHI
!    do j=GRID_JLO,GRID_JHI
!       do i=GRID_ILO,GRID_IHI
!         if((gr_meshMe .eq. MASTER_PE) .and. (abs(solnData(DUST_VAR,i,j,k)) .gt. .1)) write(*,*) i,j,k,solnData(DUST_VAR,i,j,k)
!       enddo
!    enddo
!  enddo
!   i=NGUARD+1 + 5
!   j=NGUARD+1 + 5
!   k=NGUARD+1 + 5
!   write(*,*) 'Div 1 cell=',solnData(DUST_VAR,i,j,k)
!   write(*,*) 'Uvel=',i,j,k,facexdata(VELC_FACE_VAR,i,j,k),facexdata(VELC_FACE_VAR,i+1,j,k)
!   write(*,*) 'Vvel=',i,j,k,faceydata(VELC_FACE_VAR,i,j,k),faceydata(VELC_FACE_VAR,i,j+1,k)
!   write(*,*) 'Wvel=',i,j,k,facezdata(VELC_FACE_VAR,i,j,k),facezdata(VELC_FACE_VAR,i,j,k+1),facezdata(VELC_FACE_VAR,i,j,k+2)
! endif
 

 !write(*,*) 'Mype=',gr_meshMe,gr_axisMe(IAXIS),gr_axisMe(JAXIS),gr_axisMe(KAXIS),gr_axisNumProcs(IAXIS:KAXIS)
 !call Driver_abortFlash("Done .")


#endif

  solnData(PRES_VAR,:,:,:) = 0.0
  solnData(DELP_VAR,:,:,:) = 0.0
  solnData(DUST_VAR,:,:,:) = 0.0
  solnData(TVIS_VAR,:,:,:) = 0.0

  facexData(RHDS_FACE_VAR,:,:,:) = 0.0
  faceyData(RHDS_FACE_VAR,:,:,:) = 0.0
  facezData(RHDS_FACE_VAR,:,:,:) = 0.0


  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
  call Grid_releaseBlkPtr(blockID,facezData,FACEZ)


  return

end subroutine Simulation_initBlock
