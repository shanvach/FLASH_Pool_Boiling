!! source/physics/ImBound/ImBoundMain/LevelSet
!!
!! NAME
!!
!! ib_lset(blockCount,blockList,dt)
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
!! Authors - Elizabeth Gregorio, Akash Dhruv

subroutine ib_lset(blockCount,blockList,dt)

#include "Flash.h"
#include "constants.h"
#include "SolidMechanics.h"

  ! Modules Used
  use SolidMechanics_data, only: sm_bodyInfo,sm_meshMe
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

  include "Flash_mpi.h"

  ! IO variables
  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList
  real,    INTENT(IN) :: dt

  ! Internal Variables
  integer :: numPart, e, ptelem,  nel, p
  integer, allocatable, dimension(:) :: max_ptelem

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
  real, allocatable, dimension(:,:) :: xpos,ypos
  real, allocatable, dimension(:) :: PA, PB, P1, P0, v1, v2
  real, allocatable, dimension(:,:) :: angl, dist
  real :: u,dot,det
  integer :: nelm=2,ibd ! Dimension for the points, 2 for (x,y) in 2-D
  integer :: countit
  real    :: miny, maxy, mratio, nratio, xit

  allocate(max_ptelem(gr_sbNumBodies))

  max_ptelem = 0

  do ibd=1,gr_sbNumBodies 
     if(sm_meshMe .eq. sm_BodyInfo(ibd)%BodyMaster) then
        max_ptelem(ibd) = sm_bodyInfo(ibd)%nnp-1
     end if
  end do
 
  do ibd=1,gr_sbNumBodies
     call MPI_BCAST(max_ptelem(ibd), 1, FLASH_INTEGER, sm_BodyInfo(ibd)%BodyMaster, MPI_COMM_WORLD, ierr)
  end do

  allocate(xpos(maxval(max_ptelem),gr_sbNumBodies),ypos(maxval(max_ptelem),gr_sbNumBodies))
  allocate(angl(maxval(max_ptelem),gr_sbNumBodies),dist(maxval(max_ptelem),gr_sbNumBodies))
  allocate(PA(nelm),PB(nelm),P1(nelm),P0(nelm),v1(nelm),v2(nelm))

  xpos = 0.0
  ypos = 0.0
  
  angl = 0.0
  dist = 0.0

  PA = 0.0
  PB = 0.0
  P0 = 0.0
  P1 = 0.0

  v1 = 0.0
  v2 = 0.0

  ! MPI procedure to transfer body info from local to all procs

  do ibd=1,gr_sbNumBodies
     if(sm_meshMe .eq. sm_BodyInfo(ibd)%BodyMaster) then
        xpos(1:max_ptelem(ibd),ibd) = sm_bodyInfo(ibd)%xB(2:max_ptelem(ibd)+1) + sm_bodyInfo(ibd)%qn(sm_bodyInfo(ibd)%ID(1,2:max_ptelem(ibd)+1))
        ypos(1:max_ptelem(ibd),ibd) = sm_bodyInfo(ibd)%yB(2:max_ptelem(ibd)+1) + sm_bodyInfo(ibd)%qn(sm_bodyInfo(ibd)%ID(2,2:max_ptelem(ibd)+1))
     end if
  end do

  do ibd=1,gr_sbNumBodies
        call MPI_BCAST(xpos(:,ibd), max_ptelem(ibd), FLASH_REAL, sm_BodyInfo(ibd)%BodyMaster, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(ypos(:,ibd), max_ptelem(ibd), FLASH_REAL, sm_BodyInfo(ibd)%BodyMaster, MPI_COMM_WORLD, ierr)
  end do

  ! End MPI procedure

  ! Loop through all the grid points
  do lb = 1,blockCount
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

           dist = 0.0

           ! x and y coordinates for the current grid cell
           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS)

           ycell = coord(JAXIS) - bsize(JAXIS)/2.0 +   &
                   real(j - NGUARD - 1)*del(JAXIS) +   &
                   0.5*del(JAXIS)

           zcell  = 0.0 

           do ibd=1,gr_sbNumBodies ! Loop through bodies

               countit = 0         ! Counter to check no. of intersections with 
                                   ! the body

           do p_i=1,max_ptelem(ibd) ! p_i is short for panel_index

             !if (mph_meshMe .eq. 0) print*,"Starting Lagrangian point loop"

             ! End points for the line segment of the IB
             ! PA is on the left and PB is on the right
               PA = (/xpos(p_i,ibd), ypos(p_i,ibd)/)

             !if (mph_meshMe .eq. 0) print*,"PA = (",PA,")"      

               if (p_i .eq. max_ptelem(ibd)) then
                  PB = (/xpos(1,ibd), ypos(1,ibd)/)
               else
                  PB = (/xpos(p_i+1,ibd), ypos(p_i+1,ibd)/)
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

               dist(p_i,ibd) = sqrt(v1(1)**2 + v1(2)**2)
 
               ! Find if the horizontal ray on right-side intersects with body
               miny = min(PA(2),PB(2))
               maxy = max(PA(2),PB(2))

               if(ycell .gt. miny .and. ycell .lt. maxy) then

                 ! Method #1 use ratios to divide the current panel using
                 ! y intersection and find x

                 !mratio = PA(2) - ycell
                 !nratio = ycell - PB(2)
                 !xit = (mratio*PB(1) + nratio*PA(1))/(mratio + nratio)


                 ! Method #2 use the equation of line instead

                 mratio = (PB(2)-PA(2))/(PB(1)-PA(1))          
                 xit = PA(1) + (ycell - PA(2))/mratio

                 ! Check to make sure that the intersection is on the right

                 if(xit .ge. xcell) countit = countit + 1

               end if

           end do
       
          ! Get minimum absoulte distance
          mvd = minval(dist(:,ibd))
      
          ! Construct level set - if intersections are positive then the point
          ! lies outside (-), if odd then the point lies inside (+)

          ! For first body explicitly satisfy level set, and then compare with
          ! existing level set for successive bodies

          if(ibd .eq. 1) then

               if(mod(countit,2) == 0) then
                  solnData(LMDA_VAR,i,j,k) = -mvd
               else
                  solnData(LMDA_VAR,i,j,k) = mvd
               end if
               
          else

               if(mod(countit,2) == 0) then
                  solnData(LMDA_VAR,i,j,k) = max(solnData(LMDA_VAR,i,j,k),-mvd)
               else
                  solnData(LMDA_VAR,i,j,k) = max(solnData(LMDA_VAR,i,j,k),mvd)
               endif

          end if

          end do !End body loop
         end do
        end do

        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  end do

!--------------------------------------------------------------------
    gcMask = .FALSE.
    gcMask(LMDA_VAR) = .TRUE.

    call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=ACTIVE_BLKS)
!--------------------------------------------------------------------

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
  deallocate(max_ptelem)

end subroutine ib_lset
