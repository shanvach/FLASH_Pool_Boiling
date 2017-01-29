!!****if* source/physics/IncompNS/IncompNSMain/vardens/ins_convectVelout
!!
!!
!! NAME
!!
!!  ins_convectVelout
!!
!!
!! SYNOPSIS
!!
!!  ins_convectVelout(integer(IN) :: blockCount,
!!                    integer(IN) :: blockList(blockCount)
!!                    real  (OUT) :: convvel(LOW:HIGH,MDIM))
!!
!!
!! DESCRIPTION
!!
!! Computes convective velocities out of OUTFLOW_INS domain boundaries. This routine is part 
!! of a global mass balance strategy. FOR THE MOMENT OUTFLOW_INS is defined only on HIGH x,y or z. 
!!
!! ARGUMENTS
!!
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  convvel    - Velocities out of domain for OUTFLOW_INS Boundary conditions.
!!
!!***

!#define LIQUID_PHASE_AVG 
#define TWO_PHASE_SEP

#ifdef TWO_PHASE_SEP
subroutine ins_convectVelout( blockCount, blockList, convvel1, convvel2)
#else
subroutine ins_convectVelout( blockCount, blockList, convvel)
#endif

#include "Flash.h"

  use Grid_interface, only : Grid_getDeltas,         &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_solvePoisson, Grid_getBlkBoundBox, Grid_getBlkCenterCoords

  use Grid_data, only : gr_domainBC,gr_meshComm, gr_globalDomain

#if NDIM == 2
  use IncompNS_data, only : uvel_x,vvel_x,wvel_x,uvel_y,vvel_y,wvel_y
#elif NDIM == 3  
  use IncompNS_data, only : uvel_x,vvel_x,wvel_x,uvel_y,vvel_y,wvel_y, &
                            uvel_z,vvel_z,wvel_z
#endif

  implicit none

#include "constants.h"
  include "Flash_mpi.h"


  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList 
#ifdef TWO_PHASE_SEP
  real,    INTENT(OUT) :: convvel1(LOW:HIGH,MDIM), convvel2(LOW:HIGH,MDIM)
#else
  real,    INTENT(OUT) :: convvel(LOW:HIGH,MDIM)
#endif
  !! -----------------------------------------------------

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
            
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData
  real, pointer, dimension(:,:,:,:) :: solnData

  integer :: lb,blockID,ierr,i,j,k

  integer :: faces(2,MDIM),onBoundary(2,MDIM)

  integer :: nxc,nyc,nzc

  real :: convveli(LOW:HIGH,MDIM),del(MDIM),dx,dy,dz,dxdy,dydz,dxdz,Lx,Ly,Lz
 
  real :: factorarea

  real :: convveli1(LOW:HIGH,MDIM), convveli2(LOW:HIGH,MDIM)
  real :: Ly1, Ly2, factorarea1, factorarea2
!=============================================================================

  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1
#if NDIM == 2
  nzc = -1
#elif NDIM == 3
  nzc = NZB + NGUARD + 1
#endif

#ifdef TWO_PHASE_SEP
! The free surface is assumed to be y=0 at the outlet

  convveli1 = 0.
  convveli2 = 0.
  convvel1  = 0.
  convvel2  = 0.

  ! Detect if the problem is an outflow problem to proceed with flow computation.
  ! ONLY for HIGH x,y,z.
  if (any(gr_domainBC(HIGH,1:NDIM) .eq. OUTFLOW_INS)) then    

  ! Compute Mass flow from boundaries:
  do lb = 1,blockCount
      blockID = blockList(lb)
      ! Get blocks dx, dy ,dz:
      call Grid_getDeltas(blockID,del)
      dx = del(IAXIS)
      dy = del(JAXIS)
      dxdy = dx*dy

      Lx =(gr_globalDomain(HIGH,IAXIS)-gr_globalDomain(LOW,IAXIS))
      Ly =(gr_globalDomain(HIGH,JAXIS)-gr_globalDomain(LOW,JAXIS))

     ! Shizhao modified the Ly for two phase flow outlet condition
     ! Jul 13, 2015
      Ly1 = gr_globalDomain(HIGH,JAXIS)
      Ly2 = -gr_globalDomain(LOW,JAXIS)

#if NDIM == 2
      dz = 1.
      Lz = 1.
#elif NDIM == 3
      dz = del(KAXIS)
      Lz = (gr_globalDomain(HIGH,KAXIS)-gr_globalDomain(LOW,KAXIS))
#endif
      dydz = dy*dz      
      dxdz = dx*dz
      ! Get Blocks internal limits indexes:
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 
      ! Get blocks BCs:
      call Grid_getBlkBC(blockID,faces,onBoundary)
 
      ! AXIS:
      call Grid_getBlkPtr(blockID,solnData,CENTER)
      call Grid_getBlkPtr(blockID,facexData,FACEX)
      call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
      call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

      ! High X Boundary:
      if ((faces(HIGH,IAXIS) .eq. OUTFLOW_INS)) then
!print*,"Conv A",blockID,blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS),NGUARD,nxc
        ! redistribute velocities to guardcells:
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              !Unxc-1
              uvel_x(NGUARD,j,k,HIGH,blockID) = facexData(VELC_FACE_VAR,nxc-1,j,k)  !3=19
              !Unxc
              uvel_x(NGUARD+1,j,k,HIGH,blockID) = facexData(VELC_FACE_VAR,nxc,j,k)  !4=20
           enddo
         enddo

        factorarea1 = dy*dz/(Ly1*Lz)
        factorarea2 = dy*dz/(Ly2*Lz)
        ! Compute the Mean Convective velocity across X High Boundary:
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           if(solnData(DFUN_VAR,blkLimits(HIGH,IAXIS),j,k) > 0.0d0 ) then ! air phase 
              convveli1(HIGH,IAXIS) = convveli1(HIGH,IAXIS) + & 
                                     facexData(VELC_FACE_VAR,nxc,j,k)*factorarea1
            else
              convveli2(HIGH,IAXIS) = convveli2(HIGH,IAXIS) + & 
                                     facexData(VELC_FACE_VAR,nxc,j,k)*factorarea2
            endif
            enddo
         enddo

        ! redistribute velocities for V and W:
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
              !Vnxc-1
              vvel_x(NGUARD-1,j,k,HIGH,blockID) = faceyData(VELC_FACE_VAR,nxc-1,j,k)
              !Vnxc
              vvel_x(NGUARD,j,k,HIGH,blockID)   = faceyData(VELC_FACE_VAR,nxc,j,k)
           enddo
        enddo

#if NDIM == 3
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)+1
           do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              !Wnxc-1
              wvel_x(NGUARD-1,j,k,HIGH,blockID) =  facezData(VELC_FACE_VAR,nxc-1,j,k)
              !Wnxc
              wvel_x(NGUARD,j,k,HIGH,blockID) =  facezData(VELC_FACE_VAR,nxc,j,k)
           enddo
        enddo
#endif

      endif

      ! High Y Boundary:
      if ((faces(HIGH,JAXIS) .eq. OUTFLOW_INS)) then
        ! redistribute velocities to guardcells:
!print*,"Conv B",blockID
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              !Vnyc-1
              !faceyData(VELC_FACE_VAR,i,nyc+1,k) = faceyData(VELC_FACE_VAR,i,nyc-1,k)
              vvel_y(NGUARD,i,k,HIGH,blockID) = faceyData(VELC_FACE_VAR,i,nyc-1,k)
              !Vnyc
              !faceyData(VELC_FACE_VAR,i,nyc+2,k) = faceyData(VELC_FACE_VAR,i,nyc,k) 
              vvel_y(NGUARD+1,i,k,HIGH,blockID) = faceyData(VELC_FACE_VAR,i,nyc,k)
            enddo
         enddo


        factorarea = dx*dz/(Lx*Lz)
        ! Compute the Mean Convective velocity across Y High Boundary:
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              
              convveli(HIGH,JAXIS) = convveli(HIGH,JAXIS) + & 
                                     faceyData(VELC_FACE_VAR,i,nyc,k)*factorarea
            enddo
         enddo
        

        ! redistribute velocities for U and W:
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1
              !Unyc-1
              !facexData(VELC_FACE_VAR,i,nyc+1,k) = facexData(VELC_FACE_VAR,i,nyc-1,k)
              uvel_y(NGUARD-1,i,k,HIGH,blockID) = facexData(VELC_FACE_VAR,i,nyc-1,k)
              !Unyc
              uvel_y(NGUARD,i,k,HIGH,blockID) = facexData(VELC_FACE_VAR,i,nyc,k)
           enddo
        enddo

#if NDIM == 3
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)+1
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              !Wnyc-1
              !facezData(VELC_FACE_VAR,i,nyc+1,k) = facezData(VELC_FACE_VAR,i,nyc-1,k)
              wvel_y(NGUARD-1,i,k,HIGH,blockID) = facezData(VELC_FACE_VAR,i,nyc-1,k)
              !Wnyc
              wvel_y(NGUARD,i,k,HIGH,blockID) = facezData(VELC_FACE_VAR,i,nyc,k)
           enddo
        enddo
#endif


      endif


#if NDIM == 3

      ! High Z Boundary:
      if ((faces(HIGH,KAXIS) .eq. OUTFLOW_INS)) then
        ! redistribute velocities to guardcells:
!print*,"Conv C",blockID
        do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              !Wnxc-1
              !facezData(VELC_FACE_VAR,i,j,nzc+1) = facezData(VELC_FACE_VAR,i,j,nzc-1)
              wvel_z(NGUARD,i,j,HIGH,blockID) = facezData(VELC_FACE_VAR,i,j,nzc-1)
              !Wnxc
              !facezData(VELC_FACE_VAR,i,j,nzc+2) = facezData(VELC_FACE_VAR,i,j,nzc) 
              wvel_z(NGUARD+1,i,j,HIGH,blockID) = facezData(VELC_FACE_VAR,i,j,nzc)
            enddo
         enddo


        factorarea = dx*dy/(Lx*Ly)
        ! Compute the Mean Convective velocity across Z High Boundary:
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              
              convveli(HIGH,KAXIS) = convveli(HIGH,KAXIS) + & 
                                     facezData(VELC_FACE_VAR,i,j,nzc)*factorarea
            enddo
         enddo
        

        ! redistribute velocities for U and V:
        do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1
              !Unxc-1
              !facexData(VELC_FACE_VAR,i,j,nzc+1) = facexData(VELC_FACE_VAR,i,j,nzc-1)
              uvel_z(NGUARD-1,i,j,HIGH,blockID) = facexData(VELC_FACE_VAR,i,j,nzc-1)
              !Unxc
              uvel_z(NGUARD,i,j,HIGH,blockID) =  facexData(VELC_FACE_VAR,i,j,nzc)
           enddo
        enddo

        do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              !Vnxc-1
              !faceyData(VELC_FACE_VAR,i,j,nzc+1) = faceyData(VELC_FACE_VAR,i,j,nzc-1)
              vvel_z(NGUARD-1,i,j,HIGH,blockID) = faceyData(VELC_FACE_VAR,i,j,nzc-1)
              !Vnxc
              vvel_z(NGUARD,i,j,HIGH,blockID) = faceyData(VELC_FACE_VAR,i,j,nzc)
           enddo
        enddo

      endif
#endif

      call Grid_releaseBlkPtr(blockID,solnData,CENTER)
      call Grid_releaseBlkPtr(blockID,facexData,FACEX)
      call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
      call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif


   enddo

   ! Gather total inflow volume flow ratio:
   call MPI_Allreduce(convveli1,convvel1,(HIGH-LOW+1)*MDIM,FLASH_REAL,    &
                      FLASH_SUM, gr_meshComm, ierr)
   call MPI_Allreduce(convveli2,convvel2,(HIGH-LOW+1)*MDIM,FLASH_REAL,    &
                      FLASH_SUM, gr_meshComm, ierr)


 endif ! Test if there is an OUTFLOW or NEUMANN INS BC

#else  /* TWO_PHASE_SEP*/

  convveli = 0.
  convvel  = 0.

#ifdef LIQUID_PHASE_AVG
! Just use the liquid phase to compute the convective velocity at the outlet
! This is based on the assumptions that 
     !  (1) the free surface is at y = 0 at the outlet
     !  (2) the air phase is identified by y > 0
     !  (3) the air phase has little effect on the free surface
     !  (4) the convective velocity at the outlet is 0 when y>0
! More attention should be paid when using this modification
! Jul 13, 2015
! Shizhao

  ! Detect if the problem is an outflow problem to proceed with flow computation.
  ! ONLY for HIGH x,y,z.
  if (any(gr_domainBC(HIGH,1:NDIM) .eq. OUTFLOW_INS)) then    

  ! Compute Mass flow from boundaries:
  do lb = 1,blockCount
      blockID = blockList(lb)
      ! Get blocks dx, dy ,dz:
      call Grid_getDeltas(blockID,del)
      dx = del(IAXIS)
      dy = del(JAXIS)
      dxdy = dx*dy

      Lx =(gr_globalDomain(HIGH,IAXIS)-gr_globalDomain(LOW,IAXIS))
      Ly =(gr_globalDomain(HIGH,JAXIS)-gr_globalDomain(LOW,JAXIS))

     ! Shizhao modified the Ly for two phase flow outlet condition
     ! Jul 13, 2015
      Ly = -gr_globalDomain(LOW,JAXIS)

#if NDIM == 2
      dz = 1.
      Lz = 1.
#elif NDIM == 3
      dz = del(KAXIS)
      Lz = (gr_globalDomain(HIGH,KAXIS)-gr_globalDomain(LOW,KAXIS))
#endif
      dydz = dy*dz      
      dxdz = dx*dz
      ! Get Blocks internal limits indexes:
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 
      ! Get blocks BCs:
      call Grid_getBlkBC(blockID,faces,onBoundary)
 
      ! AXIS:
      call Grid_getBlkPtr(blockID,solnData,CENTER)
      call Grid_getBlkPtr(blockID,facexData,FACEX)
      call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
      call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

      ! High X Boundary:
      if ((faces(HIGH,IAXIS) .eq. OUTFLOW_INS)) then
!print*,"Conv A",blockID,blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS),NGUARD,nxc
        ! redistribute velocities to guardcells:
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           if(solnData(DFUN_VAR,blkLimits(HIGH,IAXIS),j,k) > 0.0d0 ) then ! air phase 
              !Unxc-1
              uvel_x(NGUARD,j,k,HIGH,blockID) = 0.0d0 
              !Unxc
              uvel_x(NGUARD+1,j,k,HIGH,blockID) = 0.0d0
           else
              !Unxc-1
              uvel_x(NGUARD,j,k,HIGH,blockID) = facexData(VELC_FACE_VAR,nxc-1,j,k)  !3=19
              !Unxc
              uvel_x(NGUARD+1,j,k,HIGH,blockID) = facexData(VELC_FACE_VAR,nxc,j,k)  !4=20
           endif
           enddo
         enddo

        factorarea = dy*dz/(Ly*Lz)
        ! Compute the Mean Convective velocity across X High Boundary:
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           if(solnData(DFUN_VAR,blkLimits(HIGH,IAXIS),j,k) <= 0.0d0 ) then ! liquid phase 
              convveli(HIGH,IAXIS) = convveli(HIGH,IAXIS) + & 
                                     facexData(VELC_FACE_VAR,nxc,j,k)*factorarea
            endif
            enddo
         enddo

        ! redistribute velocities for V and W:
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
           if(solnData(DFUN_VAR,blkLimits(HIGH,IAXIS),j,k) > 0.0d0 ) then ! air phase 
              !Vnxc-1
              vvel_x(NGUARD-1,j,k,HIGH,blockID) = 0.0d0 
              !Vnxc
              vvel_x(NGUARD,j,k,HIGH,blockID)   = 0.0d0
           else
              !Vnxc-1
              vvel_x(NGUARD-1,j,k,HIGH,blockID) = faceyData(VELC_FACE_VAR,nxc-1,j,k)
              !Vnxc
              vvel_x(NGUARD,j,k,HIGH,blockID)   = faceyData(VELC_FACE_VAR,nxc,j,k)
           endif
           enddo
        enddo

#if NDIM == 3
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)+1
           do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              !Wnxc-1
              wvel_x(NGUARD-1,j,k,HIGH,blockID) =  facezData(VELC_FACE_VAR,nxc-1,j,k)
              !Wnxc
              wvel_x(NGUARD,j,k,HIGH,blockID) =  facezData(VELC_FACE_VAR,nxc,j,k)

           enddo
        enddo
#endif

      endif

      ! High Y Boundary:
      if ((faces(HIGH,JAXIS) .eq. OUTFLOW_INS)) then
        ! redistribute velocities to guardcells:
!print*,"Conv B",blockID
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              !Vnyc-1
              !faceyData(VELC_FACE_VAR,i,nyc+1,k) = faceyData(VELC_FACE_VAR,i,nyc-1,k)
              vvel_y(NGUARD,i,k,HIGH,blockID) = faceyData(VELC_FACE_VAR,i,nyc-1,k)
              !Vnyc
              !faceyData(VELC_FACE_VAR,i,nyc+2,k) = faceyData(VELC_FACE_VAR,i,nyc,k) 
              vvel_y(NGUARD+1,i,k,HIGH,blockID) = faceyData(VELC_FACE_VAR,i,nyc,k)
            enddo
         enddo


        factorarea = dx*dz/(Lx*Lz)
        ! Compute the Mean Convective velocity across Y High Boundary:
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              
              convveli(HIGH,JAXIS) = convveli(HIGH,JAXIS) + & 
                                     faceyData(VELC_FACE_VAR,i,nyc,k)*factorarea
            enddo
         enddo
        

        ! redistribute velocities for U and W:
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1
              !Unyc-1
              !facexData(VELC_FACE_VAR,i,nyc+1,k) = facexData(VELC_FACE_VAR,i,nyc-1,k)
              uvel_y(NGUARD-1,i,k,HIGH,blockID) = facexData(VELC_FACE_VAR,i,nyc-1,k)
              !Unyc
              uvel_y(NGUARD,i,k,HIGH,blockID) = facexData(VELC_FACE_VAR,i,nyc,k)
           enddo
        enddo

#if NDIM == 3
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)+1
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              !Wnyc-1
              !facezData(VELC_FACE_VAR,i,nyc+1,k) = facezData(VELC_FACE_VAR,i,nyc-1,k)
              wvel_y(NGUARD-1,i,k,HIGH,blockID) = facezData(VELC_FACE_VAR,i,nyc-1,k)
              !Wnyc
              wvel_y(NGUARD,i,k,HIGH,blockID) = facezData(VELC_FACE_VAR,i,nyc,k)
           enddo
        enddo
#endif


      endif


#if NDIM == 3

      ! High Z Boundary:
      if ((faces(HIGH,KAXIS) .eq. OUTFLOW_INS)) then
        ! redistribute velocities to guardcells:
!print*,"Conv C",blockID
        do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              !Wnxc-1
              !facezData(VELC_FACE_VAR,i,j,nzc+1) = facezData(VELC_FACE_VAR,i,j,nzc-1)
              wvel_z(NGUARD,i,j,HIGH,blockID) = facezData(VELC_FACE_VAR,i,j,nzc-1)
              !Wnxc
              !facezData(VELC_FACE_VAR,i,j,nzc+2) = facezData(VELC_FACE_VAR,i,j,nzc) 
              wvel_z(NGUARD+1,i,j,HIGH,blockID) = facezData(VELC_FACE_VAR,i,j,nzc)
            enddo
         enddo


        factorarea = dx*dy/(Lx*Ly)
        ! Compute the Mean Convective velocity across Z High Boundary:
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              
              convveli(HIGH,KAXIS) = convveli(HIGH,KAXIS) + & 
                                     facezData(VELC_FACE_VAR,i,j,nzc)*factorarea
            enddo
         enddo
        

        ! redistribute velocities for U and V:
        do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1
              !Unxc-1
              !facexData(VELC_FACE_VAR,i,j,nzc+1) = facexData(VELC_FACE_VAR,i,j,nzc-1)
              uvel_z(NGUARD-1,i,j,HIGH,blockID) = facexData(VELC_FACE_VAR,i,j,nzc-1)
              !Unxc
              uvel_z(NGUARD,i,j,HIGH,blockID) =  facexData(VELC_FACE_VAR,i,j,nzc)
           enddo
        enddo

        do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              !Vnxc-1
              !faceyData(VELC_FACE_VAR,i,j,nzc+1) = faceyData(VELC_FACE_VAR,i,j,nzc-1)
              vvel_z(NGUARD-1,i,j,HIGH,blockID) = faceyData(VELC_FACE_VAR,i,j,nzc-1)
              !Vnxc
              vvel_z(NGUARD,i,j,HIGH,blockID) = faceyData(VELC_FACE_VAR,i,j,nzc)
           enddo
        enddo

      endif
#endif

      call Grid_releaseBlkPtr(blockID,solnData,CENTER)
      call Grid_releaseBlkPtr(blockID,facexData,FACEX)
      call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
      call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif


   enddo

   ! Gather total inflow volume flow ratio:
   call MPI_Allreduce(convveli,convvel,(HIGH-LOW+1)*MDIM,FLASH_REAL,    &
                      FLASH_SUM, gr_meshComm, ierr)


 endif ! Test if there is an OUTFLOW or NEUMANN INS BC

#else  /* LIQUID_PHASE_AVG*/

!================ 
! Single phase average
    
  ! Detect if the problem is an outflow problem to proceed with flow computation.
  ! ONLY for HIGH x,y,z.
  if (any(gr_domainBC(HIGH,1:NDIM) .eq. OUTFLOW_INS)) then    

  ! Compute Mass flow from boundaries:
  do lb = 1,blockCount
      blockID = blockList(lb)
      ! Get blocks dx, dy ,dz:
      call Grid_getDeltas(blockID,del)
      dx = del(IAXIS)
      dy = del(JAXIS)
      dxdy = dx*dy

      Lx =(gr_globalDomain(HIGH,IAXIS)-gr_globalDomain(LOW,IAXIS))
      Ly =(gr_globalDomain(HIGH,JAXIS)-gr_globalDomain(LOW,JAXIS))

#if NDIM == 2
      dz = 1.
      Lz = 1.
#elif NDIM == 3
      dz = del(KAXIS)
      Lz = (gr_globalDomain(HIGH,KAXIS)-gr_globalDomain(LOW,KAXIS))
#endif
      dydz = dy*dz      
      dxdz = dx*dz
      ! Get Blocks internal limits indexes:
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 
      ! Get blocks BCs:
      call Grid_getBlkBC(blockID,faces,onBoundary)
 
      ! AXIS:
      call Grid_getBlkPtr(blockID,facexData,FACEX)
      call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
      call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

      ! High X Boundary:
      if ((faces(HIGH,IAXIS) .eq. OUTFLOW_INS)) then
!print*,"Conv A",blockID,blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS),NGUARD,nxc
        ! redistribute velocities to guardcells:
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
   
              !Unxc-1
              uvel_x(NGUARD,j,k,HIGH,blockID) = facexData(VELC_FACE_VAR,nxc-1,j,k)  !3=19
              !Unxc
              uvel_x(NGUARD+1,j,k,HIGH,blockID) = facexData(VELC_FACE_VAR,nxc,j,k)  !4=20
        
            enddo
         enddo

        factorarea = dy*dz/(Ly*Lz)
        ! Compute the Mean Convective velocity across X High Boundary:
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              
              convveli(HIGH,IAXIS) = convveli(HIGH,IAXIS) + & 
                                     facexData(VELC_FACE_VAR,nxc,j,k)*factorarea
            enddo
         enddo

        ! redistribute velocities for V and W:
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
              !Vnxc-1
              vvel_x(NGUARD-1,j,k,HIGH,blockID) = faceyData(VELC_FACE_VAR,nxc-1,j,k)
              !Vnxc
              vvel_x(NGUARD,j,k,HIGH,blockID)   = faceyData(VELC_FACE_VAR,nxc,j,k)
           enddo
        enddo

#if NDIM == 3
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)+1
           do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              !Wnxc-1
              wvel_x(NGUARD-1,j,k,HIGH,blockID) =  facezData(VELC_FACE_VAR,nxc-1,j,k)
              !Wnxc
              wvel_x(NGUARD,j,k,HIGH,blockID) =  facezData(VELC_FACE_VAR,nxc,j,k)

           enddo
        enddo
#endif

      endif

      ! High Y Boundary:
      if ((faces(HIGH,JAXIS) .eq. OUTFLOW_INS)) then
        ! redistribute velocities to guardcells:
!print*,"Conv B",blockID
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              !Vnyc-1
              !faceyData(VELC_FACE_VAR,i,nyc+1,k) = faceyData(VELC_FACE_VAR,i,nyc-1,k)
              vvel_y(NGUARD,i,k,HIGH,blockID) = faceyData(VELC_FACE_VAR,i,nyc-1,k)
              !Vnyc
              !faceyData(VELC_FACE_VAR,i,nyc+2,k) = faceyData(VELC_FACE_VAR,i,nyc,k) 
              vvel_y(NGUARD+1,i,k,HIGH,blockID) = faceyData(VELC_FACE_VAR,i,nyc,k)
            enddo
         enddo


        factorarea = dx*dz/(Lx*Lz)
        ! Compute the Mean Convective velocity across Y High Boundary:
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              
              convveli(HIGH,JAXIS) = convveli(HIGH,JAXIS) + & 
                                     faceyData(VELC_FACE_VAR,i,nyc,k)*factorarea
            enddo
         enddo
        

        ! redistribute velocities for U and W:
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1
              !Unyc-1
              !facexData(VELC_FACE_VAR,i,nyc+1,k) = facexData(VELC_FACE_VAR,i,nyc-1,k)
              uvel_y(NGUARD-1,i,k,HIGH,blockID) = facexData(VELC_FACE_VAR,i,nyc-1,k)
              !Unyc
              uvel_y(NGUARD,i,k,HIGH,blockID) = facexData(VELC_FACE_VAR,i,nyc,k)
           enddo
        enddo

#if NDIM == 3
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)+1
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              !Wnyc-1
              !facezData(VELC_FACE_VAR,i,nyc+1,k) = facezData(VELC_FACE_VAR,i,nyc-1,k)
              wvel_y(NGUARD-1,i,k,HIGH,blockID) = facezData(VELC_FACE_VAR,i,nyc-1,k)
              !Wnyc
              wvel_y(NGUARD,i,k,HIGH,blockID) = facezData(VELC_FACE_VAR,i,nyc,k)
           enddo
        enddo
#endif


      endif


#if NDIM == 3

      ! High Z Boundary:
      if ((faces(HIGH,KAXIS) .eq. OUTFLOW_INS)) then
        ! redistribute velocities to guardcells:
!print*,"Conv C",blockID
        do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              !Wnxc-1
              !facezData(VELC_FACE_VAR,i,j,nzc+1) = facezData(VELC_FACE_VAR,i,j,nzc-1)
              wvel_z(NGUARD,i,j,HIGH,blockID) = facezData(VELC_FACE_VAR,i,j,nzc-1)
              !Wnxc
              !facezData(VELC_FACE_VAR,i,j,nzc+2) = facezData(VELC_FACE_VAR,i,j,nzc) 
              wvel_z(NGUARD+1,i,j,HIGH,blockID) = facezData(VELC_FACE_VAR,i,j,nzc)
            enddo
         enddo


        factorarea = dx*dy/(Lx*Ly)
        ! Compute the Mean Convective velocity across Z High Boundary:
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              
              convveli(HIGH,KAXIS) = convveli(HIGH,KAXIS) + & 
                                     facezData(VELC_FACE_VAR,i,j,nzc)*factorarea
            enddo
         enddo
        

        ! redistribute velocities for U and V:
        do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1
              !Unxc-1
              !facexData(VELC_FACE_VAR,i,j,nzc+1) = facexData(VELC_FACE_VAR,i,j,nzc-1)
              uvel_z(NGUARD-1,i,j,HIGH,blockID) = facexData(VELC_FACE_VAR,i,j,nzc-1)
              !Unxc
              uvel_z(NGUARD,i,j,HIGH,blockID) =  facexData(VELC_FACE_VAR,i,j,nzc)
           enddo
        enddo

        do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              !Vnxc-1
              !faceyData(VELC_FACE_VAR,i,j,nzc+1) = faceyData(VELC_FACE_VAR,i,j,nzc-1)
              vvel_z(NGUARD-1,i,j,HIGH,blockID) = faceyData(VELC_FACE_VAR,i,j,nzc-1)
              !Vnxc
              vvel_z(NGUARD,i,j,HIGH,blockID) = faceyData(VELC_FACE_VAR,i,j,nzc)
           enddo
        enddo

      endif
#endif

      call Grid_releaseBlkPtr(blockID,facexData,FACEX)
      call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
      call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif


   enddo

   ! Gather total inflow volume flow ratio:
   call MPI_Allreduce(convveli,convvel,(HIGH-LOW+1)*MDIM,FLASH_REAL,    &
                      FLASH_SUM, gr_meshComm, ierr)


 endif ! Test if there is an OUTFLOW or NEUMANN INS BC

#endif  /* LIQUID_PHASE_AVG */

#endif  /* TWO_PHASE_SEP*/
 return
 end subroutine ins_convectVelout
