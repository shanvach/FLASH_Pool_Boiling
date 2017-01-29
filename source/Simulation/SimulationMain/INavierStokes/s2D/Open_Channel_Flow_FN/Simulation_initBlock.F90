!!****if* source/Simulation/SimulationMain/INavierStokes/2D/IB_Cyl_parallelIBVP_HYPRE_VD_halfDiamBel/Simulation_initBlock
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
                              sim_gCell

  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkBoundBox,    &
                             Grid_getBlkCenterCoords

  use Driver_data, ONLY : dr_simTime

!  use ins_interface, ONLY : f, fprime

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  integer, intent(in) :: blockID
  !!$ ---------------------------------
 
  integer :: i, j, k
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) ::  blIndSize,blIndSizeGC
  integer :: iii,jjj,kkk
  integer :: je,ke
  real, dimension(MDIM)  :: coord,bsize,del
  real ::  boundBox(2,MDIM)
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData

  real :: xcell,xedge,ycell,yedge,yycell,xxcell,etha_x,etha_y,temp,temp2 ! etha_y uses center coordinates for x. etha_x uses center coordinate for y.
  real :: x_blas = 428.67
  real :: etha_coeff = 13.7477    
  !----------------------------------------------------------------------

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%% f and fprime functions for blasius solution ***
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!CONTAINS
!REAL FUNCTION f (dummy_ethat)

!REAL, INTENT(IN) :: dummy_etha
!REAL :: temp
!temp = 0.0008788*dummy_etha**4 -0.02434*dummy_etha**3 +0.246*dummy_etha**2 &
!     -0.07853*dummy_etha +0.0108
!f =temp
!END FUNCTION f

!REAL FUNCTION fprime (sec_dummy_ehtha)

!REAL , INTENT(IN) :: sec_dummy_ehtha
!REAL :: temp
!temp   = -0.0000006457*sec_dummy_ehtha**8 +0.00002999*sec_dummy_ehtha**7 &
!         -0.0005358*sec_dummy_ehtha**6 +0.004655*sec_dummy_ehtha**5 &
!         -0.01963*sec_dummy_ehtha**4 + 0.03181*sec_dummy_ehtha**3 &
!         -0.02734*sec_dummy_ehtha**2 +0.3404*sec_dummy_ehtha - 0.00006472
!fprime = temp
!END FUNCTION fprime
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! write(*,*)'CALLING INITIALCONDITIONS@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@' 
  !if (myPE .eq. MASTER_PE) write(*,*) 'InitBlockTime =',dr_simTime

  ! Get nxb, nyb and nxb:
  !call Grid_getBlkIndexSize(blockId,blIndSize,blIndSizeGC)

  !nxb = blIndSize(1)
  !nyb = blIndSize(2)
  !nzb = blIndSize(3)

  ! Get Coord and Bsize for the block:
  ! Bounding box:
  call Grid_getBlkBoundBox(blockId,boundBox)
  bsize(:) = boundBox(2,:) - boundBox(1,:)

  call Grid_getBlkCenterCoords(blockId,coord)

  ! Get blocks dx, dy ,dz:
  call Grid_getDeltas(blockID,del)

  ! Point to Blocks centered variables:
  call Grid_getBlkPtr(blockID,solnData,CENTER)

  ! Point to Blocks face variables: 
  call Grid_getBlkPtr(blockID,facexData,FACEX)
  call Grid_getBlkPtr(blockID,faceyData,FACEY)

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)


  !- kpd - Initialize the distance function in the 1st quadrant 
  do k=1,blkLimitsGC(HIGH,KAXIS)
     do j=1,blkLimitsGC(HIGH,JAXIS)
        do i=1,blkLimitsGC(HIGH,IAXIS)

           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS)

           ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS)  +  &
                   0.5*del(JAXIS)

          !if (ycell .LE. 0.0) then
          !   solnData(DFUN_VAR,i,j,k) = 0.0 - ycell
          !else
          !   solnData(DFUN_VAR,i,j,k) = -1.0 * ycell
          !end if
          !solnData(DFUN_VAR,i,j,k) = ycell - 2*(2.7182**(-0.1*xcell))*sin(2*xcell)
         solnData(DFUN_VAR,i,j,k) = ycell !-sin(xcell)
        enddo
     enddo
  enddo



  ! set values for u,v velocities and pressure
  solnData(PRES_VAR,:,:,:) = 0.0
  solnData(DELP_VAR,:,:,:) = 0.0
  solnData(DUST_VAR,:,:,:) = 0.0
  solnData(TVIS_VAR,:,:,:) = 0.0

  solnData(CURV_VAR,:,:,:) = 0.0
  solnData(SIGP_VAR,:,:,:) = 0.0
  solnData(VISC_VAR,:,:,:) = 0.0
  solnData(PFUN_VAR,:,:,:) = 0.0
  !solnData(DFUN_VAR,:,:,:) = 0.0

 ! facexData(VELC_FACE_VAR,:,:,:) = 1.0
!#################################################
!#### Setting X and Y-face values ################
!#################################################

  do k=1,blkLimitsGC(HIGH,KAXIS)
        do i=1,blkLimitsGC(HIGH,IAXIS)
            do j=1,blkLimitsGC(HIGH,JAXIS)
           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS)
           xxcell = xcell - 0.5*del(IAXIS)
           ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS)  +  &
                   0.5*del(JAXIS)
           yycell = ycell - 0.5*del(JAXIS)
          
          if (yycell .gt. 0.0) then
       
                facexData(VELC_FACE_VAR,i,j,k) = 0.0
                faceyData(VELC_FACE_VAR,i,j,k) = 0.0
          else
                etha_x = (etha_coeff*abs(ycell))/sqrt(xxcell + x_blas)
                etha_y = (etha_coeff*abs(yycell))/sqrt(xcell + x_blas)

     !   temp2=  -0.0000006457*etha_x**8 +0.00002999*etha_x**7 &
     !           -0.0005358*etha_x**6 +0.004655*etha_x**5 &
     !           -0.01963*etha_x**4 +0.03181*etha_x**3 &
     !           -0.02734*etha_x**2 +0.3404*etha_x -0.00006472



          if ( etha_x .lt. 8.0) then
        
        facexData(VELC_FACE_VAR,i,j,k) =1.0*(-0.0000006457*etha_x**8 +0.00002999*etha_x**7 &
                                             -0.0005358*etha_x**6 +0.004655*etha_x**5 &
                                             -0.01963*etha_x**4 + 0.03181*etha_x**3 &
                                             -0.02734*etha_x**2 +0.3404*etha_x - 0.00006472)

        temp =  -0.0000006457*etha_y**8 +0.00002999*etha_y**7 &
                -0.0005358*etha_y**6 +0.004655*etha_y**5 &
                -0.01963*etha_y**4 + 0.03181*etha_y**3 &
                -0.02734*etha_y**2 +0.3404*etha_y - 0.00006472       
 
        temp2=  0.0008788*etha_y**4 -0.02434*etha_y**3 +0.246*etha_y**2 &
                -0.07853*etha_y +0.0108

        faceyData(VELC_FACE_VAR,i,j,k) = -1.0*(0.5*sqrt(0.0056497/(xcell+x_blas))*(etha_y*temp-temp2))
        else
                temp =0.0
                temp2=0.0

         facexData(VELC_FACE_VAR,i,j,k)=1.0*(-0.0000006457*(8.0**8) +0.00002999*(8.0**7) &
                                             -0.0005358*(8.0**6) +0.004655*(8.0**5) &
                                             -0.01963*(8.0**4) + 0.03181*(8.0**3) &
                                             -0.02734*(8.0**2) +0.3404*8.0 - 0.00006472)

        temp = (-0.0000006457*(8.0**8) +0.00002999*(8.0**7) &
         -0.0005358*(8.0**6) +0.004655*(8.0**5) &
         -0.01963*(8.0**4) + 0.03181*(8.0**3) &
         -0.02734*(8.0**2) +0.3404*8.0 - 0.00006472)

        temp2=  0.0008788*(8.0**4) -0.02434*(8.0**3)+0.246*(8.0**2) &
                -0.07853*8.0 +0.0108

        faceyData(VELC_FACE_VAR,i,j,k) =-1.0*(0.5*sqrt(0.0056497/(xcell+x_blas))*(8.0*temp-temp2)) 
       endif
        endif          
        enddo
     enddo
  enddo

!##########################################################


  facexData(RHDS_FACE_VAR,:,:,:) = 0.0
  faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

  facexData(SIGM_FACE_VAR,:,:,:) = 0.0
  faceyData(SIGM_FACE_VAR,:,:,:) = 0.0
  facexData(RH1F_FACE_VAR,:,:,:) = 0.0
  faceyData(RH1F_FACE_VAR,:,:,:) = 0.0
  facexData(RH2F_FACE_VAR,:,:,:) = 0.0
  faceyData(RH2F_FACE_VAR,:,:,:) = 0.0


!!$  ! Point to blocks center and face vars:
!!$  call Grid_getBlkPtr(blockID,solnData,CENTER)
!!$  call Grid_getBlkPtr(blockID,facexData,FACEX)
!!$  call Grid_getBlkPtr(blockID,faceyData,FACEY)
!!$
!!$
!!$  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
!!$     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
!!$     do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
!!$
!!$
!!$     if (ISNAN(facexData(VELC_FACE_VAR,i,j,k))) then
!!$       write(*,*) 'facexData block=',blockID
!!$       write(*,*) 'i,j,k=',i,j,k,' is a NAN.',facexData(VELC_FACE_VAR,i,j,k)
!!$     endif
!!$
!!$     if (ISNAN(faceyData(VELC_FACE_VAR,i,j,k))) then
!!$       write(*,*) 'faceyData block=',blockID
!!$       write(*,*) 'i,j,k=',i,j,k,' is a NAN.',faceyData(VELC_FACE_VAR,i,j,k)
!!$     endif
!!$
!!$
!!$     enddo
!!$     enddo
!!$  enddo

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)



  return

111    format (i4,3x,i4)
112    format (3(3x,e12.4))

end subroutine Simulation_initBlock
