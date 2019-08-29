!! source/physics/Multiphase/MultiphaseMain/INS/ImBound
!!
!! NAME
!!
!! mph_iblset(blockCount,blockList)
!!
!! SYNOPSIS
!! 
!!  
!! VARIABLES
!!
!!
!! DESCRIPTION
!! 
!! Subroutine to find the distance function lambda for 
!! the immersed boundary (IB).
!! 
!! Author - Elizabeth Gregorio 

subroutine mph_iblset(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

#include "Flash.h"
#include "constants.h"
#include "SolidMechanics.h"

  ! Modules Used
  use SolidMechanics_data, only: sm_bodyInfo
  use sm_element_interface, only: sm_el02_mapParticles, sm_el10_mapParticles, &
                                  sm_el01_mapParticles

  use gr_sbData, ONLY : gr_sbBodyInfo,gr_sbNumBodies
  use Driver_interface, ONLY : Driver_abortFlash

  use Grid_interface, ONLY : Grid_getListOfBlocks,   &
                             Grid_getDeltas,         &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_fillGuardCells,    &
                             Grid_getBlkBoundBox,Grid_getBlkCenterCoords

  use Multiphase_data, only: mph_rho1,mph_rho2,mph_sten,mph_crmx,mph_crmn, &
                             mph_vis1,mph_vis2,mph_lsit, mph_inls, mph_meshMe

  use ib_interface, ONLY : ib_stencils
  use ImBound_data, only : ib_stencil
!  use gr_sbData, only: xb, yb

  implicit none

  ! IO variables
  integer, intent(in) :: sweepOrder
  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList
  real,    INTENT(IN) :: timeEndAdv,dt,dtOld

  ! Internal Variables
  integer :: numPart, e, ptelem, max_ptelem, nel, p, max_ptm1
  !integer, allocatable, dimension(:) :: loc_num

  ! Arugments List
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox

  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS), isAttached

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  integer :: lb,ii,jj,kk,ierr,i,j,k,dir,blockID,p_i,pp_i

  real :: mva, mvd

  real bsize(MDIM),coord(MDIM)

  real del(MDIM),xcell,ycell,zcell

  integer :: listofBlocks(MAXBLOCKS)
  integer :: intval

  ! For the algorithm
  real, allocatable, dimension(:) :: xpos,ypos
  real, allocatable, dimension(:) :: PA, PB, P1, P0, v1, v2
  real, allocatable, dimension(:) :: angl, dist
  !integer, allocatable, dimension(:):: ind
  real :: u,dot,det
  integer :: nelm=2,ibd ! Dimension for the points, 2 for (x,y) in 2-D

  !integer :: max_ptelem=4 ! for box should be 4

  do ibd=1,gr_sbNumBodies
      max_ptelem = sm_bodyInfo(ibd)%nnp
      !max_ptelem = nel
      max_ptm1 = max_ptelem-1
      MPI_Bcast()
  end do

  allocate(xpos(max_ptelem),ypos(max_ptelem))
  allocate(angl(max_ptelem),dist(max_ptelem))
  allocate(PA(nelm),PB(nelm),P1(nelm),P0(nelm),v1(nelm),v2(nelm))
  !allocate(ind(nelm))

  do ibd=1,gr_sbNumBodies
      xpos = sm_bodyInfo(ibd)%xB!(1:nel)
      ypos = sm_bodyInfo(ibd)%yB!(1:nel)
      MPI_Bcast()
      MPI_Bcast()
  end do

  !xpos = (/3.0,  -3.0, -3.0,  3.0/)
  !ypos = (/-3.0, -3.0, -3.25, -3.25/)

  ! Loop through all the grid points
  do lb = 1,blockCount
        if (mph_meshMe .eq. 0) print*,"Starting blockCount loop"
        ! Loop through all the blocks
        blockID = blockList(lb)

        ! Get the block Id and boundaries
        call Grid_getBlkBoundBox(blockId,boundBox)

        ! Size of the block (far side - near side)
        bsize(:) = boundBox(2,:) - boundBox(1,:)

        ! Get co-ordinates of the block:
        call Grid_getBlkCenterCoords(blockId,coord)

        ! Get the delta x and y for the block:
        call Grid_getDeltas(blockID,del)
        
        ! Get Blocks internal limits indexes:
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

        ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
        call Grid_getBlkPtr(blockID,facezData,FACEZ)
        
        k = 1
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
         do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           if (mph_meshMe .eq. 0) print*,"Starting Eulerian grid loop"    
           angl = 0.0
           dist = 0.0

           ! x and y coordinates for the current grid cell
           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS)

           ycell = coord(JAXIS) - bsize(JAXIS)/2.0 +   &
                   real(j - NGUARD - 1)*del(JAXIS) +   &
                   0.5*del(JAXIS)

           zcell  = 0.0 

           do p_i=1,max_ptelem ! p_i is short for panel_index

             !if (mph_meshMe .eq. 0) print*,"Starting Lagrangian point loop"

             ! End points for the line segment of the IB
             ! PA is on the left and PB is on the right
               PA = (/xpos(p_i), ypos(p_i)/)

             !if (mph_meshMe .eq. 0) print*,"PA = (",PA,")"      

               if (p_i .eq. max_ptelem) then
                  PB = (/xpos(1), ypos(1)/)
               else
                  PB = (/xpos(p_i+1), ypos(p_i+1)/)
               end if
             !if (mph_meshMe .eq. 0) print*,"PB = (",PB,")"      

!             if (mph_meshMe .eq. 0) print*,"PB assigned, p_i = ",p_i,", max_ptelem = ",max_ptelem

             ! Grid cell point
               P1 = (/xcell, ycell/)

             ! Drop a normal from P1 to the line made by connecting PA PB (not the
             ! line segment)
               u = ((P1(1)-PA(1))*(PB(1)-PA(1)) + (P1(2)-PA(2))*(PB(2)-PA(2))) / &
                   (((PB(1)-PA(1))**2)+((PB(2)-PA(2))**2)) 
        
             ! Re-assign u if the normal hits the line segment to the left of PA or
             ! the right of PB
               if (u .lt. 0) then
                  u = 0.0
               else if (u .gt. 1) then
                  u = 1.0
               end if 

             ! Find the point on the line segment with the shortest distance to P1
             ! (If the normal hits the line outside the line segment it is
             !  reassigned to hit the closer endpoint.)
               P0(1) = PA(1) + (PB(1) - PA(1))*u
               P0(2) = PA(2) + (PB(2) - PA(2))*u

             ! Determine the quadrent and angle for the "normal"
             ! (If to the left or right of the line segment the vector with the 
             !  shortest distance to the line segment will not be perpendicular)
             
               if (abs(P0(1)-PA(1)) .lt. 1e-13 .and. abs(P0(2)-PA(2)) .lt. 1e-13) then
                  v1 = (/(P1(1) - P0(1)),(P1(2) - P0(2))/)
                  v2 = (/(P0(1) - PB(1)),(P0(2) - PB(2))/)
               else
                  v1 = (/(P1(1) - P0(1)),(P1(2) - P0(2))/)
                  v2 = (/(PA(1) - P0(1)),(PA(2) - P0(2))/)
               end if

               dot =   v1(1)*v2(1) + v1(2)*v2(2)
               det = -(v1(1)*v2(2) - v1(2)*v2(1))
  
               angl(p_i) = atan2(det, dot)
               dist(p_i) = sqrt(v1(1)**2 + v1(2)**2)

           end do


           !if ( dist(1) .lt. dist(2) ) then
               !mvd = dist!(1)
               !mva = angl!(1)
           !else
           !    mvd = dist(2)
           !    mva = angl(2)
           !end if

           !ind = minloc(dist) ! Finds the index of the minimum distance
           !mvd = dist(ind(1)) ! Gets the minimum distance value
           !mva = angl(ind(1)) ! Gets the angle for the minimum distance
         
           do pp_i=1,max_ptelem 
           do p_i=1,max_ptelem-1
              if (dist(p_i) > dist(p_i+1)) then

                       dist(p_i) = dist(p_i) + dist(p_i+1)
                       dist(p_i+1) = dist(p_i) - dist(p_i+1)
                       dist(p_i) = dist(p_i) - dist(p_i+1)

                       angl(p_i) = angl(p_i) + angl(p_i+1)
                       angl(p_i+1) = angl(p_i) - angl(p_i+1)
                       angl(p_i) = angl(p_i) - angl(p_i+1)

                end if
           end do
           end do

           mvd = dist(1)
           mva = angl(1)

  !        if (mph_meshMe .eq. 0) print*,"Minimum Values, mvd = ",mvd,", mva = ",mva
 
           if (mva .eq. 0.0) then
               solnData(LMDA_VAR,i,j,k) = 0.0
           else
               solnData(LMDA_VAR,i,j,k) = -mvd*sign(1.0,mva)
           end if

         end do
        end do

        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  end do

  print*,"Loops ended"

!--------------------------------------------------------------------
    gcMask = .FALSE.
    gcMask(LMDA_VAR) = .TRUE.

    call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=ACTIVE_BLKS)
!--------------------------------------------------------------------

  if (mph_meshMe .eq. 0) print*,"masks done"

  do lb=1,blockCount

        ! Loop through all the blocks
        blockID = blockList(lb)

        ! Get the block Id and boundaries
        call Grid_getBlkBoundBox(blockId,boundBox)

        ! Size of the block (far side - near side)
        bsize(:) = boundBox(2,:) - boundBox(1,:)

        ! Get co-ordinates of the block:
        call Grid_getBlkCenterCoords(blockId,coord)

        ! Get the delta x and y for the block:
        call Grid_getDeltas(blockID,del)
        
        ! Get Blocks internal limits indexes:
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

        ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
        call Grid_getBlkPtr(blockID,facezData,FACEZ)
        
        k = 1
        do j=2,blkLimitsGC(HIGH,JAXIS)-1
            do i=2,blkLimitsGC(HIGH,IAXIS)-1

               solnData(NMLX_VAR,i,j,k) = -((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))/&
                                      sqrt(((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))**2)

               solnData(NMLY_VAR,i,j,k) = -((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(IAXIS))/&
                                      sqrt(((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))**2)


            end do
        end do

        !k = 1
        do j=2,blkLimitsGC(HIGH,JAXIS)-1
            do i=2,blkLimitsGC(HIGH,IAXIS)-1
               solnData(TNGY_VAR,i,j,k) = ((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))/&
                                     sqrt(((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))**2+&
                                          ((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))**2)

               solnData(TNGX_VAR,i,j,k) = -((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(IAXIS))/&
                                      sqrt(((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))**2)

            end do
        end do        

        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  end do

  deallocate(xpos,ypos)
  deallocate(angl,dist)
  deallocate(PA,PB,P1,P0,v1,v2)
  !deallocate(ind)

end subroutine mph_iblset
