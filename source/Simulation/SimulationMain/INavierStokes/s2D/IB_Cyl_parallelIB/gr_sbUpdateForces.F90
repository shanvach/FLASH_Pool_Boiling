!!****if* source/Simulation/SimulationMain/INavierStokes/2D/IB_Cyl_parallelIB/gr_sbUpdateForces
!!
!! NAME
!!  gr_sbUpdateForces
!!
!! SYNOPSIS
!!
!!  gr_sbUpdateForces()
!!
!! DESCRIPTION
!!
!! ARGUMENTS
!!
!!***


#include "constants.h"
#include "Flash.h"

subroutine gr_sbUpdateForces
  use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_meshNumProcs
  use Grid_interface, ONLY : Grid_updateSolidBodyForces, Grid_getBlkPtr, &
                             Grid_releaseBlkPtr, &
                             Grid_getDeltas, Grid_getBlkCenterCoords, &
                             Grid_getBlkPhysicalSize,Grid_getBlkIndexLimits

  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, totalPart, &
       gr_sbDebug, gr_sbParticleCount, solid_body
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use ImBound_data, only : ib_dt

  use tree, only : neigh, lrefine, lrefine_max

  implicit none
#include "Flash_mpi.h"
  type(solid_body), pointer :: bodyInfo
  real, dimension(MDIM) :: particleposn
  integer :: i, p, gettingFrom, blkID, b, recvCount
  real :: particleData(NPART_PROPS)
  real ::  Body_Cen(NDIM),Voli
  real :: MomArmi_x, MomArmi_y, Fxtot(gr_sbNumBodies), Fytot(gr_sbNumBodies), Momz(gr_sbNumBodies), &
                                Fxtoti(gr_sbNumBodies),Fytoti(gr_sbNumBodies),Momzi(gr_sbNumBodies)
  real :: Cdrag,Clift,Cmoment

  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData,facexData2,faceyData2,facezData2

  integer, save, dimension(MAXBLOCKS) :: listOfBlocks
  integer, save :: count
  integer :: lb,blockID,ierr,j,k

  real :: FxTotL, FyTotL, MomzL, Fx, Fy, FxProc, FyProc, MomZProc
  real :: del(MDIM),coord(MDIM),bsize(MDIM),dx,dy,dxdy,xcell,ycell

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  integer :: low_x_low_y, low_x_high_y, high_x_low_y, high_x_high_y

  call Timers_start("update_forces")

  call Grid_getListOfBlocks(LEAF, listOfBlocks, count)
  do lb=1,count

     blockID = listOfBlocks(lb)

     ! Get face data (velocities):
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

     facexData(FORC_FACE_VAR,:,:,:) = 0.
     faceyData(FORC_FACE_VAR,:,:,:) = 0.

#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
     facezData(FORC_FACE_VAR,:,:,:) = 0.
#endif

     ! Release face data (velocities):
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

  enddo


  ! Loop over all bodies
  do b = 1, gr_sbNumBodies
     bodyInfo => gr_sbBodyInfo(b)

     if (bodyInfo % myPE == bodyInfo % bodyMaster) then

        ! Initialize Particle Lagrangian Forces:
        bodyInfo % particles(FUL_PART_PROP,:) = 0.
        bodyInfo % particles(FVL_PART_PROP,:) = 0.
#if NDIM == 3
        bodyInfo % particles(FWL_PART_PROP,:) = 0.
#endif

        do i = 1, totalPart
           if (int(bodyInfo % particles(PROC_PART_PROP,i)) == bodyInfo % bodyMaster) then
              particleposn(IAXIS) = (bodyInfo % particles(POSX_PART_PROP,i))
              if(NDIM >1) then
                 particleposn(JAXIS) = (bodyInfo % particles(POSY_PART_PROP,i))
              endif
              if(NDIM >2) then
                 particleposn(KAXIS) = (bodyInfo % particles(POSZ_PART_PROP,i))
              endif
              blkID = int(bodyInfo % particles(BLK_PART_PROP,i))
              particleData = bodyInfo % particles(1:NPART_PROPS,i)
              call Grid_updateSolidBodyForces(blkID, particleData)
              !call Grid_updateSolidBodyForces(blkID, b, i, particleposn)
              bodyInfo % particles(1:NPART_PROPS,i) = particleData   

              if (i .eq. 1) then
                 ! The Block center coord
                 call Grid_getBlkCenterCoords(blkID,coord)
     
                 ! The Block Physical size
                 call Grid_getBlkPhysicalSize(blkID,bsize)

                 write(*,*) i,bodyInfo%particles(POSX_PART_PROP,i),bodyInfo%particles(POSY_PART_PROP,i)
                 write(*,*) blkID,coord(1)-bsize(1)/2.,coord(1)+bsize(1)/2.,coord(2)-bsize(2)/2.,coord(2)+bsize(2)/2. 
              endif

           end if
        end do

     else
        
        gettingFrom = gr_sbParticleCount(b)
        if (gettingFrom > 0) then

        ! Initialize Particle Lagrangian Forces:
        bodyInfo % particles(FUL_PART_PROP,:) = 0.
        bodyInfo % particles(FVL_PART_PROP,:) = 0.
#if NDIM == 3
        bodyInfo % particles(FWL_PART_PROP,:) = 0.
#endif


           recvCount = gettingFrom
           
           do p = 1, recvCount
              blkID = int(bodyInfo % particles(BLK_PART_PROP,p))
              particleposn(IAXIS) = bodyInfo % particles(POSX_PART_PROP,p)
              if(NDIM >1) then
                 particleposn(JAXIS) = bodyInfo % particles(POSY_PART_PROP,p)
              endif
              if(NDIM >2) then
                 particleposn(KAXIS) = bodyInfo % particles(POSZ_PART_PROP,p)
              endif
              particleData =  bodyInfo % particles(1:NPART_PROPS,p)
              call Grid_updateSolidBodyForces(blkID, particleData)
              !updating the forces

              bodyInfo % particles(1:NPART_PROPS,p) = particleData 
              !call Grid_updateSolidBodyForces(blkID, b, p, particleposn)
              !write(*,'(a, i6, a,3f8.2,a,i6)') "myID", gr_meshMe, "receiving position", & 
              !RecvBuf(POSX_PART_PROP:POSZ_PART_PROP,p), "from", gr_sbBodyInfo(b) % bodyMaster 
           enddo
        
          !deallocate(bodyInfo % particles)
        end if
     end if

  end do  ! End Loop over Bodies


  ! Do inverse guardcell filling:
  do lb=1,count

     blockID = listOfBlocks(lb)

     ! Get face data (velocities):
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

!     if (any(facexData(FORC_FACE_VAR,:,:,:) .ne. 0.)) &
!     write(*,*) blockID
!     if (blockID .eq. 202) &
!     write(*,*) facexData(FORC_FACE_VAR,:,:,:)

     facexData(FORB_FACE_VAR,:,:,:) = facexData(FORC_FACE_VAR,:,:,:)

     faceyData(FORB_FACE_VAR,:,:,:) = faceyData(FORC_FACE_VAR,:,:,:)

     ! Release face data (velocities):
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

  enddo

  do lb=1,count

     blockID = listOfBlocks(lb)

     ! Get face data (velocities):
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)


     if ((lrefine(blockID) .eq. lrefine_max)) then     

     if (neigh(1,1,blockID) .gt. 0) then ! Lower x

        ! Get face data (velocities):
        call Grid_getBlkPtr(neigh(1,1,blockID),facexData2,FACEX)
        call Grid_getBlkPtr(neigh(1,1,blockID),faceyData2,FACEY)

        ! force in xdir
        facexData(FORC_FACE_VAR,NGUARD+1:NGUARD+2,NGUARD+1:NGUARD+NYB,:) = &
        facexData(FORC_FACE_VAR,NGUARD+1:NGUARD+2,NGUARD+1:NGUARD+NYB,:) + &
        facexData2(FORB_FACE_VAR,NGUARD+NXB+1:NGUARD+NXB+2,NGUARD+1:NGUARD+NYB,:) ! GuardCell vals

        ! force in ydir
        faceyData(FORC_FACE_VAR,NGUARD:NGUARD+1,NGUARD+1:NGUARD+NYB+1,:) = &
        faceyData(FORC_FACE_VAR,NGUARD:NGUARD+1,NGUARD+1:NGUARD+NYB+1,:) + &
        faceyData2(FORB_FACE_VAR,NGUARD+NXB:NGUARD+NXB+1,NGUARD+1:NGUARD+NYB+1,:) ! GuardCell vals
    

        ! Release face data (velocities):
        call Grid_releaseBlkPtr(neigh(1,1,blockID),facexData2,FACEX)
        call Grid_releaseBlkPtr(neigh(1,1,blockID),faceyData2,FACEY)
     endif

     if (neigh(1,2,blockID) .gt. 0) then ! Higher x

        ! Get face data (velocities):
        call Grid_getBlkPtr(neigh(1,2,blockID),facexData2,FACEX)
        call Grid_getBlkPtr(neigh(1,2,blockID),faceyData2,FACEY)

        ! force in xdir
        facexData(FORC_FACE_VAR,NGUARD+NXB:NGUARD+NXB+1,NGUARD+1:NGUARD+NYB,:) = &
        facexData(FORC_FACE_VAR,NGUARD+NXB:NGUARD+NXB+1,NGUARD+1:NGUARD+NYB,:) + &
        facexData2(FORB_FACE_VAR,NGUARD:NGUARD+1,NGUARD+1:NGUARD+NYB,:) ! GuardCell vals

        ! force in ydir
        faceyData(FORC_FACE_VAR,NGUARD+NXB:NGUARD+NXB+1,NGUARD+1:NGUARD+NYB+1,:) = &
        faceyData(FORC_FACE_VAR,NGUARD+NXB:NGUARD+NXB+1,NGUARD+1:NGUARD+NYB+1,:) + &
        faceyData2(FORB_FACE_VAR,NGUARD:NGUARD+1,NGUARD+1:NGUARD+NYB+1,:) ! GuardCell vals
    

        ! Release face data (velocities):
        call Grid_releaseBlkPtr(neigh(1,2,blockID),facexData2,FACEX)
        call Grid_releaseBlkPtr(neigh(1,2,blockID),faceyData2,FACEY)
     endif

     if (neigh(1,3,blockID) .gt. 0) then ! Lower y

        ! Get face data (velocities):
        call Grid_getBlkPtr(neigh(1,3,blockID),facexData2,FACEX)
        call Grid_getBlkPtr(neigh(1,3,blockID),faceyData2,FACEY)

        ! force in xdir
        facexData(FORC_FACE_VAR,NGUARD+1:NGUARD+NXB+1,NGUARD:NGUARD+1,:) = &
        facexData(FORC_FACE_VAR,NGUARD+1:NGUARD+NXB+1,NGUARD:NGUARD+1,:) + &
        facexData2(FORB_FACE_VAR,NGUARD+1:NGUARD+NXB+1,NGUARD+NYB:NGUARD+NYB+1,:) ! GuardCell vals

        ! force in ydir
        faceyData(FORC_FACE_VAR,NGUARD+1:NGUARD+NXB,NGUARD+1:NGUARD+2,:) = &
        faceyData(FORC_FACE_VAR,NGUARD+1:NGUARD+NXB,NGUARD+1:NGUARD+2,:) + &
        faceyData2(FORB_FACE_VAR,NGUARD+1:NGUARD+NXB,NGUARD+NYB+1:NGUARD+NYB+2,:) ! GuardCell vals
    

        ! Release face data (velocities):
        call Grid_releaseBlkPtr(neigh(1,3,blockID),facexData2,FACEX)
        call Grid_releaseBlkPtr(neigh(1,3,blockID),faceyData2,FACEY)
     endif

     if (neigh(1,4,blockID) .gt. 0) then ! Higher y

        ! Get face data (velocities):
        call Grid_getBlkPtr(neigh(1,4,blockID),facexData2,FACEX)
        call Grid_getBlkPtr(neigh(1,4,blockID),faceyData2,FACEY)

        ! force in xdir
        facexData(FORC_FACE_VAR,NGUARD+1:NGUARD+NXB+1,NGUARD+NYB:NGUARD+NYB+1,:) = &
        facexData(FORC_FACE_VAR,NGUARD+1:NGUARD+NXB+1,NGUARD+NYB:NGUARD+NYB+1,:) + &
        facexData2(FORB_FACE_VAR,NGUARD+1:NGUARD+NXB+1,NGUARD:NGUARD+1,:) ! GuardCell vals

        ! force in ydir
        faceyData(FORC_FACE_VAR,NGUARD+1:NGUARD+NXB,NGUARD+NYB:NGUARD+NYB+1,:) = &
        faceyData(FORC_FACE_VAR,NGUARD+1:NGUARD+NXB,NGUARD+NYB:NGUARD+NYB+1,:) + &
        faceyData2(FORB_FACE_VAR,NGUARD+1:NGUARD+NXB,NGUARD:NGUARD+1,:) ! GuardCell vals
    

        ! Release face data (velocities):
        call Grid_releaseBlkPtr(neigh(1,4,blockID),facexData2,FACEX)
        call Grid_releaseBlkPtr(neigh(1,4,blockID),faceyData2,FACEY)
     endif


     ! Corners
     if (neigh(1,1,blockID) .gt. 0) then
     if (neigh(1,3,neigh(1,1,blockID)) .gt. 0) then ! Low x, low y

        low_x_low_y = neigh(1,3,neigh(1,1,blockID))

        ! Get face data (velocities):
        call Grid_getBlkPtr(low_x_low_y,facexData2,FACEX)
        call Grid_getBlkPtr(low_x_low_y,faceyData2,FACEY)

        ! force in xdir
        facexData(FORC_FACE_VAR,NGUARD+1:NGUARD+2,NGUARD:NGUARD+1,:) = &
        facexData(FORC_FACE_VAR,NGUARD+1:NGUARD+2,NGUARD:NGUARD+1,:) + &
        facexData2(FORB_FACE_VAR,NGUARD+NXB+1:NGUARD+NXB+2,NGUARD+NYB:NGUARD+NYB+1,:) ! GuardCell vals

        ! force in ydir
        faceyData(FORC_FACE_VAR,NGUARD:NGUARD+1,NGUARD+1:NGUARD+2,:) = &
        faceyData(FORC_FACE_VAR,NGUARD:NGUARD+1,NGUARD+1:NGUARD+2,:) + &
        faceyData2(FORB_FACE_VAR,NGUARD+NXB:NGUARD+NXB+1,NGUARD+NYB+1:NGUARD+NYB+2,:) ! GuardCell vals
    

        ! Release face data (velocities):
        call Grid_releaseBlkPtr(low_x_low_y,facexData2,FACEX)
        call Grid_releaseBlkPtr(low_x_low_y,faceyData2,FACEY)        

     endif


     if (neigh(1,4,neigh(1,1,blockID)) .gt. 0) then ! Low x, high y

        low_x_high_y = neigh(1,4,neigh(1,1,blockID))

        ! Get face data (velocities):
        call Grid_getBlkPtr(low_x_high_y,facexData2,FACEX)
        call Grid_getBlkPtr(low_x_high_y,faceyData2,FACEY)

        ! force in xdir
        facexData(FORC_FACE_VAR,NGUARD+1:NGUARD+2,NGUARD+NYB:NGUARD+NYB+1,:) = &
        facexData(FORC_FACE_VAR,NGUARD+1:NGUARD+2,NGUARD+NYB:NGUARD+NYB+1,:) + &
        facexData2(FORB_FACE_VAR,NGUARD+NXB+1:NGUARD+NXB+2,NGUARD:NGUARD+1,:) ! GuardCell vals

        ! force in ydir
        faceyData(FORC_FACE_VAR,NGUARD:NGUARD+1,NGUARD+NYB:NGUARD+NYB+1,:) = &
        faceyData(FORC_FACE_VAR,NGUARD:NGUARD+1,NGUARD+NYB:NGUARD+NYB+1,:) + &
        faceyData2(FORB_FACE_VAR,NGUARD+NXB:NGUARD+NXB+1,NGUARD:NGUARD+1,:) ! GuardCell vals
    

        ! Release face data (velocities):
        call Grid_releaseBlkPtr(low_x_high_y,facexData2,FACEX)
        call Grid_releaseBlkPtr(low_x_high_y,faceyData2,FACEY)        

     endif
     endif

     if (neigh(1,2,blockID) .gt. 0) then
     if (neigh(1,3,neigh(1,2,blockID)) .gt. 0) then ! High x, low y

        high_x_low_y = neigh(1,3,neigh(1,2,blockID))

        ! Get face data (velocities):
        call Grid_getBlkPtr(high_x_low_y,facexData2,FACEX)
        call Grid_getBlkPtr(high_x_low_y,faceyData2,FACEY)

        ! force in xdir
        facexData(FORC_FACE_VAR,NGUARD+NXB:NGUARD+NXB+1,NGUARD:NGUARD+1,:) = &
        facexData(FORC_FACE_VAR,NGUARD+NXB:NGUARD+NXB+1,NGUARD:NGUARD+1,:) + &
        facexData2(FORB_FACE_VAR,NGUARD:NGUARD+1,NGUARD+NYB:NGUARD+NYB+1,:) ! GuardCell vals

        ! force in ydir
        faceyData(FORC_FACE_VAR,NGUARD+NXB:NGUARD+NXB+1,NGUARD+1:NGUARD+2,:) = &
        faceyData(FORC_FACE_VAR,NGUARD+NXB:NGUARD+NXB+1,NGUARD+1:NGUARD+2,:) + &
        faceyData2(FORB_FACE_VAR,NGUARD:NGUARD+1,NGUARD+NYB+1:NGUARD+NYB+2,:) ! GuardCell vals
    

        ! Release face data (velocities):
        call Grid_releaseBlkPtr(high_x_low_y,facexData2,FACEX)
        call Grid_releaseBlkPtr(high_x_low_y,faceyData2,FACEY)        

     endif


     if (neigh(1,4,neigh(1,2,blockID)) .gt. 0) then ! High x, high y

        high_x_high_y = neigh(1,4,neigh(1,2,blockID))

        ! Get face data (velocities):
        call Grid_getBlkPtr(high_x_high_y,facexData2,FACEX)
        call Grid_getBlkPtr(high_x_high_y,faceyData2,FACEY)

        ! force in xdir
        facexData(FORC_FACE_VAR,NGUARD+NXB:NGUARD+NXB+1,NGUARD+NYB:NGUARD+NYB+1,:) = &
        facexData(FORC_FACE_VAR,NGUARD+NXB:NGUARD+NXB+1,NGUARD+NYB:NGUARD+NYB+1,:) + &
        facexData2(FORB_FACE_VAR,NGUARD:NGUARD+1,NGUARD:NGUARD+1,:) ! GuardCell vals

        ! force in ydir
        faceyData(FORC_FACE_VAR,NGUARD+NXB:NGUARD+NXB+1,NGUARD+NYB:NGUARD+NYB+1,:) = &
        faceyData(FORC_FACE_VAR,NGUARD+NXB:NGUARD+NXB+1,NGUARD+NYB:NGUARD+NYB+1,:) + &
        faceyData2(FORB_FACE_VAR,NGUARD:NGUARD+1,NGUARD:NGUARD+1,:) ! GuardCell vals
    

        ! Release face data (velocities):
        call Grid_releaseBlkPtr(high_x_high_y,facexData2,FACEX)
        call Grid_releaseBlkPtr(high_x_high_y,faceyData2,FACEY)        

     endif
     endif

     endif

     ! Release face data (velocities):
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

  enddo

  ! Do the sum of forces to Velocities
  do lb=1,count

     blockID = listOfBlocks(lb)

     ! Get face data (velocities):
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)


     facexData(VELC_FACE_VAR,:,:,:) =  facexData(VELC_FACE_VAR,:,:,:) + &
                                       ib_dt*facexData(FORC_FACE_VAR,:,:,:)
     faceyData(VELC_FACE_VAR,:,:,:) =  faceyData(VELC_FACE_VAR,:,:,:) + &
                                       ib_dt*faceyData(FORC_FACE_VAR,:,:,:)


#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
     facezData(VELC_FACE_VAR,:,:,:) =  facezData(VELC_FACE_VAR,:,:,:) + &
                                       ib_dt*facezData(FORC_FACE_VAR,:,:,:)
#endif

     ! Release face data (velocities):
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo


  ! Calculate the Lagrangian Forces on each body  
  do b = 1, gr_sbNumBodies
     bodyInfo => gr_sbBodyInfo(b)

     Fxtot(b) = 0.
     Fytot(b) = 0.
     Momz(b)  = 0.

     Fxtoti(b) = 0.
     Fytoti(b) = 0.
     Momzi(b)  = 0.

#if NDIM == 3
     `Fztot(b) = 0.; Fztoti(b) = 0.
#endif 

     ! Body Center
     Body_Cen(IAXIS) = sum(bodyInfo % xb(1:bodyInfo%NumVertices-1))/real(bodyInfo%NumVertices-1)
     Body_Cen(JAXIS) = sum(bodyInfo % yb(1:bodyInfo%NumVertices-1))/real(bodyInfo%NumVertices-1)

     if (bodyInfo % myPE == bodyInfo % bodyMaster) then
        do i = 1, totalPart
           if (int(bodyInfo % particles(PROC_PART_PROP,i)) == bodyInfo % bodyMaster) then

              MomArmi_x = bodyInfo%particles(POSY_PART_PROP,i)-Body_Cen(JAXIS)
              MomArmi_y = bodyInfo%particles(POSX_PART_PROP,i)-Body_Cen(IAXIS)

              Voli =  bodyInfo%particles(AREA_PART_PROP,i)*bodyInfo%particles(HL_PART_PROP,i)


              Fxtoti(b) =  Fxtoti(b) + bodyInfo%particles(FUL_PART_PROP,i)* Voli
              Fytoti(b) =  Fytoti(b) + bodyInfo%particles(FVL_PART_PROP,i)* Voli
              Momzi (b) =  Momzi(b) - bodyInfo%particles(FUL_PART_PROP,i)*MomArmi_x*Voli +  &
                                      bodyInfo%particles(FVL_PART_PROP,i)*MomArmi_y*Voli

#if NDIM == 3
              Fztoti(b) = Fztoti(b) + bodyInfo%particles(FWL_PART_PROP,i)* Voli
#endif 
           endif
        enddo

     else !MyPe not BodyMaster

        gettingFrom = gr_sbParticleCount(b)
        if (gettingFrom > 0) then
          recvCount = gettingFrom
           do p = 1, recvCount

              MomArmi_x = bodyInfo%particles(POSY_PART_PROP,p)-Body_Cen(JAXIS)
              MomArmi_y = bodyInfo%particles(POSX_PART_PROP,p)-Body_Cen(IAXIS)

              Voli =  bodyInfo%particles(AREA_PART_PROP,p)*bodyInfo%particles(HL_PART_PROP,p)

              Fxtoti(b) =  Fxtoti(b) + bodyInfo%particles(FUL_PART_PROP,p)* Voli
              Fytoti(b) =  Fytoti(b) + bodyInfo%particles(FVL_PART_PROP,p)* Voli
              Momzi (b) =  Momzi(b) - bodyInfo%particles(FUL_PART_PROP,p)*MomArmi_x*Voli +  &
                                      bodyInfo%particles(FVL_PART_PROP,p)*MomArmi_y*Voli


#if NDIM == 3
              Fztoti(b) = Fztoti(b) + bodyInfo%particles(FWL_PART_PROP,p)* Voli
#endif 
           enddo
        endif
     endif

     ! Sum force values among processors for each body
     ! Fx:
     call MPI_Allreduce(Fxtoti(b),Fxtot(b), CONSTANT_ONE, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)     

     ! Fy:
     call MPI_Allreduce(Fytoti(b),Fytot(b), CONSTANT_ONE, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)   

     ! Mz:
     call MPI_Allreduce(Momzi(b),Momz(b), CONSTANT_ONE, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)   


     Cdrag    = 2.* FxTot(b)
     Clift    = 2.* FyTot(b)
     Cmoment  = 2.* Momz(b)

     if (gr_meshMe .eq. MASTER_PE)  then
     write(*,*) ' '
     write(*,'(A35,3g18.12,A10,I8)') ' Total Lagr Fx , Fy & Mz=',FxTot(b), FyTot(b), Momz(b),'on body ',b
     write(*,'(A35,3g18.12)') ' Euler Drag , Lift & Moment ceoff=',Cdrag, Clift,Cmoment
     endif

  enddo


  ! Calculate the Total Eulerian Forces for all bodies:
  ! loop over the leaf blocks
  Fx = 0.; Fy=0.; MomzL=0.;
  do lb = 1,count

     blockID = listofblocks(lb)
     
     ! Get face data (velocities):
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif 
     
     ! Get dx,dy
     call Grid_getDeltas(blockID,del)
     
     ! The Block center coord
     call Grid_getBlkCenterCoords(blockId,coord)
     
     ! The Block Physical size
     call Grid_getBlkPhysicalSize(blockID,bsize)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 
     
     ! Local cell size:
     dx = del(IAXIS)
     dy = del(JAXIS)
     dxdy=dx*dy
     
     ! Calculating the total forces 

     ! X - Direction
     Fx = SUM (facexData(FORC_FACE_VAR,                 &
          blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)    &
          ))*dxdy + Fx 
     
     ! Y - Direction
     Fy = SUM (faceyData( FORC_FACE_VAR,                &
          blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)    &
          ))*dxdy + Fy      
     
     ! Calc the moment about the Z (2D) (+ve is anti clock-wise)
     do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)

!!$           yface = coord(JAXIS) - 0.5*bsize(JAXIS) + &
!!$                   real(j-NGUARD-1)*dy
   
           ycell = coord(JAXIS) - 0.5*bsize(JAXIS) + &
                   real(j-NGUARD-1)*dy + 0.5*dy

           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

!!$              xface = coord(IAXIS) - 0.5*bsize(IAXIS) + &
!!$                      real(i-NGUARD-1)*dx

              xcell = coord(IAXIS) - 0.5*bsize(IAXIS) + &
                      real(i-NGUARD-1)*dx + 0.5*dx

              MomzL = -dxdy*facexData(FORC_FACE_VAR,i,j,k)*(ycell-0.) + &
                       dxdy*faceyData(FORC_FACE_VAR,i,j,k)*(xcell-0.) + &
                       MomzL

           enddo
        enddo
     enddo

     ! Release face data (Forces):
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
     
     
  enddo
  
  FxProc = Fx
  call MPI_Allreduce(FxProc, FxTotL, 1, FLASH_REAL,&
       MPI_SUM, MPI_COMM_WORLD, ierr)
  
  FyProc = Fy  
  call MPI_Allreduce(FyProc, FyTotL, 1, FLASH_REAL,&
       MPI_SUM, MPI_COMM_WORLD, ierr)

  Momzproc = MomzL
  call MPI_Allreduce(Momzproc,MomzL, 1, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)
 
  ! Calculating the drag and lift, and Moment coefficients
  Cdrag  = 2.* FxTotL  ! Reference length, velocity and density == 1
  Clift  = 2.* FyTotL
  Cmoment= 2.* MomZL

  if (gr_meshMe .eq. MASTER_PE)  then
     write(*,*) ' '
     write(*,'(A35,3g18.12)') ' Total Euler Fx , Fy & Mz=',FxTotL, FyTotL, MomzL
     write(*,'(A35,3g18.12)') ' Euler Drag , Lift & Moment coeff=',Cdrag, Clift,Cmoment
  endif


  deallocate(gr_sbParticleCount)
  call Timers_stop("update_forces")
end subroutine gr_sbUpdateForces
