!! source/physics/ImBound/ImBoundMain/LevelSet
!!
!! NAME
!!
!! ib_lset_3D(blockCount,blockList,dt)
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
!! 3D immersed boundary with triangular nodes.
!! 
!! (This subroutine requires optimization)
!!
!! ---- Akash Dhruv

subroutine ib_lset_3D(blockCount,blockList,dt)

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

  implicit none

  include "Flash_mpi.h"

  ! IO variables
  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList
  real,    INTENT(IN) :: dt

  ! Internal Variables
  integer :: numPart, e, ptelem,  nel, p
  integer, allocatable, dimension(:) :: max_ptelem, max_wsnel

  ! Arugments List
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox

  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS), isAttached

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  integer :: lb,ii,jj,kk,ierr,i,j,k,dir,blockID,ielem

  real :: mva, mvd

  real bsize(MDIM),coord(MDIM)

  real del(MDIM),xcell,ycell,zcell

  integer :: listofBlocks(MAXBLOCKS)
  integer :: intval

  ! For the algorithm
  real, allocatable, dimension(:,:) :: xpos,ypos,zpos
  real, allocatable, dimension(:,:) :: xcenter,ycenter,zcenter
  real, allocatable, dimension(:,:,:) :: normal
  integer, allocatable, dimension(:,:,:) :: elem
  real, dimension(3) :: PA, PB, P1, P0, PC, nrm, lx, ln 
  real, dimension(3) :: min_vec, max_vec, tempvec, PP, vec, PN
  real, dimension(3) :: vecA, vecB, vecW
  real :: dotD, da, db
  real :: tempnorm, tempnorm1, tempnorm2, tempnorm3
  real, allocatable, dimension(:,:) :: dist
  real :: du,dn
  integer :: nelm=3,ibd ! Dimension for the points, 2 for (x,y) in 2-D
  integer :: countit
  integer TAIB(2),count_rateIB
  real*8  ETIB

  allocate(max_ptelem(gr_sbNumBodies))
  allocate(max_wsnel(gr_sbNumBodies))

  max_ptelem = 0
  max_wsnel = 0

  do ibd=1,gr_sbNumBodies 
     if(sm_meshMe .eq. sm_BodyInfo(ibd)%BodyMaster) then
        max_ptelem(ibd) = sm_bodyInfo(ibd)%nnp
        max_wsnel(ibd)  = sm_bodyInfo(ibd)%ws_nel
     end if
  end do

  do ibd=1,gr_sbNumBodies
     call MPI_BCAST(max_ptelem(ibd), 1, FLASH_INTEGER, sm_BodyInfo(ibd)%BodyMaster, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(max_wsnel(ibd),  1, FLASH_INTEGER, sm_BodyInfo(ibd)%BodyMaster, MPI_COMM_WORLD, ierr)
  end do

  allocate(xpos(maxval(max_ptelem),gr_sbNumBodies),ypos(maxval(max_ptelem),gr_sbNumBodies),zpos(maxval(max_ptelem),gr_sbNumBodies))
  allocate(elem(3,maxval(max_wsnel),gr_sbNumBodies))
  allocate(dist(maxval(max_wsnel),gr_sbNumBodies))
  allocate(xcenter(maxval(max_wsnel),gr_sbNumBodies),ycenter(maxval(max_wsnel),gr_sbNumBodies),zcenter(maxval(max_wsnel),gr_sbNumBodies))
  allocate(normal(3,maxval(max_wsnel),gr_sbNumBodies))

  xpos = 0.0
  ypos = 0.0
  zpos = 0.0
  elem = 0.0
  
  dist = 0.0

  PA = 0.0
  PB = 0.0
  PC = 0.0

  P0 = 0.0
  P1 = 0.0

  nrm = 0.0
  lx  = 0.0
  ln  = 0.0

  ! MPI procedure to transfer body info from local to all procs
  do ibd=1,gr_sbNumBodies
     if(sm_meshMe .eq. sm_BodyInfo(ibd)%BodyMaster) then
        xpos(1:max_ptelem(ibd),ibd) = sm_bodyInfo(ibd)%xB + sm_bodyInfo(ibd)%qn(sm_bodyInfo(ibd)%ID(1,:))
        ypos(1:max_ptelem(ibd),ibd) = sm_bodyInfo(ibd)%yB + sm_bodyInfo(ibd)%qn(sm_bodyInfo(ibd)%ID(2,:))
        zpos(1:max_ptelem(ibd),ibd) = sm_bodyInfo(ibd)%zB + sm_bodyInfo(ibd)%qn(sm_bodyInfo(ibd)%ID(3,:))
       elem(:,1:max_wsnel(ibd),ibd) = sm_bodyInfo(ibd)%ws_IEN
     end if
  end do

  do ibd=1,gr_sbNumBodies
        call MPI_BCAST(xpos(:,ibd), max_ptelem(ibd), FLASH_REAL, sm_BodyInfo(ibd)%BodyMaster, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(ypos(:,ibd), max_ptelem(ibd), FLASH_REAL, sm_BodyInfo(ibd)%BodyMaster, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(zpos(:,ibd), max_ptelem(ibd), FLASH_REAL, sm_BodyInfo(ibd)%BodyMaster, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(elem(:,:,ibd), max_wsnel(ibd)*3, FLASH_INTEGER, sm_BodyInfo(ibd)%BodyMaster, MPI_COMM_WORLD, ierr)
  end do
  ! End MPI procedure

  do ibd=1,gr_sbNumBodies
     do ielem=1,max_wsnel(ibd)

        xcenter(ielem,ibd) = (xpos(elem(1,ielem,ibd),ibd) + xpos(elem(2,ielem,ibd),ibd) + xpos(elem(3,ielem,ibd),ibd))/3
        ycenter(ielem,ibd) = (ypos(elem(1,ielem,ibd),ibd) + ypos(elem(2,ielem,ibd),ibd) + ypos(elem(3,ielem,ibd),ibd))/3
        zcenter(ielem,ibd) = (zpos(elem(1,ielem,ibd),ibd) + zpos(elem(2,ielem,ibd),ibd) + zpos(elem(3,ielem,ibd),ibd))/3

        PA = (/xpos(elem(1,ielem,ibd),ibd), ypos(elem(1,ielem,ibd),ibd), zpos(elem(1,ielem,ibd),ibd)/) 
        PB = (/xpos(elem(2,ielem,ibd),ibd), ypos(elem(2,ielem,ibd),ibd), zpos(elem(2,ielem,ibd),ibd)/) 
        PC = (/xpos(elem(3,ielem,ibd),ibd), ypos(elem(3,ielem,ibd),ibd), zpos(elem(3,ielem,ibd),ibd)/)

        call cross_product(tempvec,PB-PA,PC-PA)
        call norm(tempnorm,tempvec)

        normal(1,ielem,ibd) = tempvec(1)/tempnorm
        normal(2,ielem,ibd) = tempvec(2)/tempnorm
        normal(3,ielem,ibd) = tempvec(3)/tempnorm

     end do
  end do

  if(sm_meshMe .eq. 0) print *,"IB Level Set Reconstruction in Progess"

  CALL SYSTEM_CLOCK(TAIB(1),count_rateIB)

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
        
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
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

           zcell  = coord(KAXIS) - bsize(KAXIS)/2.0 +   &
                   real(k - NGUARD - 1)*del(KAXIS) +   &
                   0.5*del(KAXIS)

           do ibd=1,gr_sbNumBodies ! Loop through bodies

               countit = 0         ! Counter to check no. of intersections with 
                                   ! the body

           do ielem=1,max_wsnel(ibd) ! Loop through elements on each body

                P1 = (/xcell, ycell, zcell/)

                P0 = (/xcenter(ielem,ibd), ycenter(ielem,ibd), zcenter(ielem,ibd)/)
                PA = (/xpos(elem(1,ielem,ibd),ibd), ypos(elem(1,ielem,ibd),ibd), zpos(elem(1,ielem,ibd),ibd)/)
                PB = (/xpos(elem(2,ielem,ibd),ibd), ypos(elem(2,ielem,ibd),ibd), zpos(elem(2,ielem,ibd),ibd)/)
                PC = (/xpos(elem(3,ielem,ibd),ibd), ypos(elem(3,ielem,ibd),ibd), zpos(elem(3,ielem,ibd),ibd)/)

                nrm = (/normal(1,ielem,ibd), normal(2,ielem,ibd), normal(3,ielem,ibd)/)

                min_vec = (/minval((/PA(1), PB(1), PC(1)/)), minval((/PA(2), PB(2), PC(2)/)), minval((/PA(3), PB(3), PC(3)/))/)
                max_vec = (/maxval((/PA(1), PB(1), PC(1)/)), maxval((/PA(2), PB(2), PC(2)/)), maxval((/PA(3), PB(3), PC(3)/))/)

                lx = (/1.0, 0.0, 0.0/)
                ln = nrm

                vec  = PA-P1
                vecA = PB-PA;
                vecB = PC-PA;
                dotD = dot_product(vecA,vecB)**2 - dot_product(vecA,vecA)*dot_product(vecB,vecB)

                if(abs(dot_product(lx,nrm)) .lt. 1e-13) then
                   du = 0.0

                else
                   du = dot_product(vec,nrm)/dot_product(lx,nrm)

                end if
                
                PP   = P1 + du*lx
                vecW = PP-PA;                 
                da = (dot_product(vecA,vecB)*dot_product(vecW,vecB) - dot_product(vecB,vecB)*dot_product(vecW,vecA))/dotD
                db = (dot_product(vecA,vecB)*dot_product(vecW,vecA) - dot_product(vecA,vecA)*dot_product(vecW,vecB))/dotD;


                if(da .ge. 0.0 .and. da .le. 1.0 .and. db .ge. 0.0 .and. (da+db) .le. 1.0 .and. du .gt. 0.0) &
                   countit = countit + 1

                dn = dot_product(vec,nrm)/dot_product(ln,nrm)
                PN = P1 + dn*ln
                vecW = PN-PA;                 
                da = (dot_product(vecA,vecB)*dot_product(vecW,vecB) - dot_product(vecB,vecB)*dot_product(vecW,vecA))/dotD
                db = (dot_product(vecA,vecB)*dot_product(vecW,vecA) - dot_product(vecA,vecA)*dot_product(vecW,vecB))/dotD;

                if(da .ge. 0.0 .and. da .le. 1.0 .and. db .ge. 0.0 .and. (da+db) .le. 1.0) then

                        call norm(tempnorm,P1-PN)
                        dist(ielem,ibd) = tempnorm

                else

                        call norm(tempnorm1, P1-PA)
                        call norm(tempnorm2, P1-PB)
                        call norm(tempnorm3, P1-PC)
                        call norm(tempnorm,  P1-P0)

                        dist(ielem,ibd) = minval((/tempnorm, tempnorm1, tempnorm2, tempnorm3/))

                end if 

           end do
       
           ! Get minimum absoulte distance
           mvd = minval(dist(:,ibd))
      
           ! Construct level set - if intersections are positive then the point
           ! lies outside (-), if odd then the point lies inside (+)

           ! For first body explicitly satisfy level set, and then compare with
           ! existing level set for successive bodies

           if(ibd .eq. 1) then

               if(mod(countit,2) == 1) then
                  solnData(LMDA_VAR,i,j,k) = mvd
               else
                  solnData(LMDA_VAR,i,j,k) = -mvd
               end if
               
           else

               if(mod(countit,2) == 1) then
                  solnData(LMDA_VAR,i,j,k) = max(solnData(LMDA_VAR,i,j,k),mvd)
               else
                  solnData(LMDA_VAR,i,j,k) = max(solnData(LMDA_VAR,i,j,k),-mvd)
               endif

           end if

           end do !End body loop

          end do
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

        do k=2,blkLimitsGC(HIGH,KAXIS)-1
        do j=2,blkLimitsGC(HIGH,JAXIS)-1
        do i=2,blkLimitsGC(HIGH,IAXIS)-1

           solnData(NMLX_VAR,i,j,k) = -((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))/&
                                      sqrt(((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j,k+1) - solnData(LMDA_VAR,i,j,k-1))/2*del(KAXIS))**2)

           solnData(NMLY_VAR,i,j,k) = -((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))/&
                                      sqrt(((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j,k+1) - solnData(LMDA_VAR,i,j,k-1))/2*del(KAXIS))**2)

           solnData(NMLZ_VAR,i,j,k) = -((solnData(LMDA_VAR,i,j,k+1) - solnData(LMDA_VAR,i,j,k-1))/2*del(KAXIS))/&
                                      sqrt(((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j,k+1) - solnData(LMDA_VAR,i,j,k-1))/2*del(KAXIS))**2)

        end do
        end do
        end do

        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  end do

  deallocate(xpos,ypos,zpos,elem)
  deallocate(dist)
  deallocate(max_ptelem)
  deallocate(max_wsnel)
  deallocate(xcenter,ycenter,zcenter,normal)
 
  CALL SYSTEM_CLOCK(TAIB(2),count_rateIB)
  ETIB=REAL(TAIB(2)-TAIB(1),8)/count_rateIB

  if(sm_meshMe .eq. 0) print *,"IB Level Set Reconstruction Complete"
  if(sm_meshMe .eq. 0) write(*,*) 'Total Reconstruction Time =',ETIB

end subroutine ib_lset_3D

subroutine cross_product(cross,a,b)

        implicit none
        real, dimension(3), intent(out) :: cross
        real, dimension(3), intent(in) :: a, b

        cross(1) = a(2) * b(3) - a(3) * b(2)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)        

end subroutine cross_product

subroutine norm(mag,vec)

        implicit none
        real, intent(out) :: mag
        real, dimension(3), intent(in)  :: vec
        
        mag = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)

end subroutine norm
