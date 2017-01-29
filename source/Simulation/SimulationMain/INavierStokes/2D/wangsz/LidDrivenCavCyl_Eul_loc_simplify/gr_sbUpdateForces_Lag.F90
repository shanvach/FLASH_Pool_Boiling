!!****if* source/physics/ImBound/ImBoundMain/LagForce/Extras/gr_sbUpdateForces
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

#define DO_MARKER_FORCING /* If commented, no IB forcing from Markers */
!#define ITER_FORCING /* Iterative procedure in Kermpe and Frohlich 2005, JCP */

#include "constants.h"
#include "Flash.h"
#include "ImBound.h"

subroutine gr_sbUpdateForces
  
  use Grid_interface, ONLY : Grid_updateSolidBodyForces, Grid_getBlkPtr,     &
                             Grid_releaseBlkPtr, Grid_getListOfBlocks,       &
                             Grid_getDeltas, Grid_getBlkCenterCoords,        &
                             Grid_getBlkPhysicalSize,Grid_getBlkIndexLimits, &
                             Grid_fillGuardCells

  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, totalPart, &
                        gr_sbDebug, gr_sbParticleCount, solid_body

  use gr_ptData, ONLY :  gr_ptBlkList, gr_ptBlkCount

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use ImBound_data, only : ib_dt, ib_maxIterForcing

  use gr_ptInterface, ONLY : gr_ptMove, gr_ptSetIndices, gr_ptResetIndices

  use Particles_data, ONLY : pt_indexList, pt_indexCount,pt_maxPerProc,pt_posinitialized

  use ib_interface, only : ib_forceInsideBody

  use Driver_interface, only : Driver_getNStep

  use ins_interface, only : ins_fluxfix

#ifdef FLASH_GRID_PARAMESH
  use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_meshNumProcs, gr_maxParticlesPerProc
  use tree, only : neigh, lrefine, lrefine_max
#else
  use Grid_data, ONLY : gr_axisComm, gr_exch, gr_gridDataStruct, &
       gr_justExchangedGC,gr_domainBC, &
       gr_offset, gr_meshMe, gr_meshComm, gr_meshNumProcs
#endif


  implicit none
#include "Flash_mpi.h"
  type(solid_body), pointer :: bodyInfo
  real, dimension(MDIM) :: particleposn
  integer :: i,j,k,b,p, gettingFrom, blockID, recvCount
  real :: particleData(NPART_PROPS)
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData,facexData2,faceyData2,facezData2
  integer, save, dimension(MAXBLOCKS) :: listOfBlocks
  integer, save :: count
  integer :: lb,ierr
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: ct,localCount,globalCount,localPart
  logical :: moveDone,coords_in_blks=.true.

  real :: del(MDIM),coord(MDIM),bsize(MDIM)

  integer :: iForceIter, maxForceIter
  character(len=2) :: indForceIter
  character(len=6) :: indNStep
  integer :: NStep

  !! Forcing Test Variables
  real :: Body_Cen(MDIM),Voli
  real :: MomArmi_x, MomArmi_y, MomArmi_z
  real :: Fxtot(gr_sbNumBodies),  Fytot(gr_sbNumBodies), &
          Fxtoti(gr_sbNumBodies),Fytoti(gr_sbNumBodies)
  real :: Momz(gr_sbNumBodies),Momzi(gr_sbNumBodies)                                
  real :: Momx(gr_sbNumBodies),  Momxi(gr_sbNumBodies),  &
          Momy(gr_sbNumBodies),  Momyi(gr_sbNumBodies)
  real :: Fztot(gr_sbNumBodies), Fztoti(gr_sbNumBodies) 

  real :: dx,dy,dz,xcell,ycell,zcell,dxdydz
#ifdef FLASH_GRID_PARAMESH
  integer :: nxc,nyc,nzc
#endif

  !! Inverse gcfill variables:
#ifndef INVERSE_GCELL_FILL 

#ifdef FLASH_GRID_PARAMESH
  real :: particles2(NPART_PROPS,gr_maxParticlesPerProc)
#else
  real :: particles2(NPART_PROPS,pt_maxPerProc)
#endif

  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)

#else
#endif

  real :: FxTotL, FyTotL, FzTotL, MomxL, MomyL, MomzL, &
          Fx, Fy, Fz, FxProc, FyProc, FzProc, MomXProc, MomYProc, MomZProc


  integer, save :: countf = 0
  integer :: countRP, countVP
  character(len=28) :: fileName
  character(len=6) :: index_count,index_proc

  integer :: val
  integer :: Particles_Allbods


#ifdef FLASH_GRID_PARAMESH
  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1
  nzc = NZB + NGUARD + 1
#endif

  Particles_Allbods = 0

  !FzTotL = 0.

!  call Timers_start("update_forces")

!  call Timers_start('InitEulerForce')

#ifdef FLASH_GRID_PARAMESH
  call Grid_getListOfBlocks(LEAF, listOfBlocks, count)
  call Grid_getListOfBlocks(LEAF, gr_ptBlkList, gr_ptBlkCount)
#else
  count=CONSTANT_ONE
  listOfBlocks(CONSTANT_ONE)=CONSTANT_ONE
#endif

  Particles_Allbods = 0

  ! Set Eulerian Force Variable to zero:
  do lb=1,count

     blockID = listOfBlocks(lb)

     ! Get face data (velocities):
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

     facexData(FORC_FACE_VAR,:,:,:) = 0.
     faceyData(FORC_FACE_VAR,:,:,:) = 0.

#if NDIM == MDIM
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
     facezData(FORC_FACE_VAR,:,:,:) = 0.
#endif

     ! Release face data (velocities):
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

  enddo

!  call Timers_stop('InitEulerForce') ! Shizhao 

#ifdef DO_MARKER_FORCING
!  call Timers_start("Force_particles")

  ! Loop over all bodies
  do b = 1, gr_sbNumBodies
     bodyInfo => gr_sbBodyInfo(b)

#ifndef INVERSE_GCELL_FILL
!     call Timers_start("Create_Particles2")
     particles2(:,:) = 0.
#ifdef FLASH_GRID_PARAMESH
     globalCount = gr_maxParticlesPerProc  
     moveDone=.false.
#else
     globalCount=pt_maxPerProc 
#endif

     !! Now get all the indices into the data structure setup right
     !! for the rest of the unit
     call gr_ptSetIndices(pt_indexList, pt_indexCount)

     if (bodyInfo % myPE == bodyInfo % bodyMaster) then

!        call Timers_start('SetParticles2')

        totalPart = bodyInfo % totalPart

        ct = 0 
        do i = 1, totalPart
           if (int(bodyInfo % particles(PROC_PART_PROP,i)) == bodyInfo % bodyMaster) then
              ct= ct+1
              particles2(:,ct) = bodyInfo % particles(:,i)
           endif
        enddo
        
        localCount  = ct

!        call Timers_stop('SetParticles2')

        ! Generate virtual Particles: This call to gr_ptMove is calles with the
        ! only purpose of generating virtual particles for surface markers on the
        ! framework of the Immersed Boundary method. This might have to change if
        ! Different types of particles are used.
#ifdef FLASH_GRID_PARAMESH
        call gr_ptMove(particles2,NPART_PROPS, localCount,globalCount, moveDone)
#else
!        call Timers_start('GridMoveParticles')

!        write(*,*) 'M1, GlobalCount, LocalCount, ID:', globalCount, localCount, gr_meshMe
        pt_posinitialized=.true.
        call Grid_moveParticles(particles2,NPART_PROPS, globalCount, localCount, &
             pt_indexList, pt_indexCount, coords_in_blks)
!        write(*,*) 'M2, GlobalCount, LocalCount, ID:', globalCount, localCount, gr_meshMe
!        call Timers_stop('GridMoveParticles')
#endif
     else

!        call Timers_start('SetParticles2')
        localCount = gr_sbParticleCount(b)
  
        if (localCount .gt. 0) &
        particles2(:,1:localCount) = bodyInfo%particles(:,1:localCount)

!        call Timers_stop('SetParticles2')

        ! Generate virtual Particles:
#ifdef FLASH_GRID_PARAMESH     
        call gr_ptMove(particles2,NPART_PROPS, localCount,globalCount, moveDone)
#else
!        call Timers_start('GridMoveParticles')
!        write(*,*) 'S1, GlobalCount, LocalCount, ID:', globalCount, localCount, gr_meshMe
        pt_posinitialized=.true.
        call Grid_moveParticles(particles2,NPART_PROPS, globalCount, localCount, &
             pt_indexList, pt_indexCount, coords_in_blks)
!        write(*,*) 'S2, GlobalCount, LocalCount, ID:', globalCount, localCount, gr_meshMe
!        call Timers_stop('GridMoveParticles')
#endif
     endif
!     call Timers_stop("Create_Particles2")

     call gr_ptResetIndices(pt_indexList, pt_indexCount)     

#endif  /* #ifndef INVERSE_GCELL_FILL */


     recvCount = 0
     if (bodyInfo % myPE == bodyInfo % bodyMaster) then

#ifndef INVERSE_GCELL_FILL

        ! Find number of particles in the Master Processor:
        localPart = 0
        do j=1,size(bodyInfo%particlesPerProc,DIM=2)
           if(bodyInfo%particlesPerProc(1,j) .eq. bodyInfo%bodyMaster) then
              localPart = bodyInfo%particlesPerProc(2,j)
           endif
        enddo

        particles2(FUL_PART_PROP,:) = 0.
        particles2(FVL_PART_PROP,:) = 0.
#if NDIM == MDIM
        particles2(FWL_PART_PROP,:) = 0.
#endif
        recvCount = localCount

        do i = 1, localCount

           if (int(particles2(PROC_PART_PROP,i)) == bodyInfo % bodyMaster) then
#ifdef FLASH_GRID_PARAMESH
              blockID = int(particles2(BLK_PART_PROP,i))
#else
              blockID=CONSTANT_ONE
#endif
              particleData = particles2(1:NPART_PROPS,i)


#else /* #ifndef INVERSE_GCELL_FILL */

#endif /* #ifndef INVERSE_GCELL_FILL */

              !call Timers_start("Grid_updateSolidBodyForces")      
              call Grid_updateSolidBodyForces(b,int(particleData(GLOB_PART_PROP)),blockID, particleData)
              !call Timers_stop("Grid_updateSolidBodyForces")

#ifndef INVERSE_GCELL_FILL
              if (i .le. localPart) &
              bodyInfo % particles(1:NPART_PROPS,i) = particleData   
#else
#endif

           end if
        end do

#ifdef INVERSE_GCELL_FILL
#endif

     else ! if (bodyInfo % myPE == bodyInfo % bodyMaster) then 
        
        gettingFrom = gr_sbParticleCount(b)

#ifndef INVERSE_GCELL_FILL 

        if (localCount > 0) then
           particles2(FUL_PART_PROP,:) = 0.
           particles2(FVL_PART_PROP,:) = 0.
#if NDIM == MDIM
           particles2(FWL_PART_PROP,:) = 0.
#endif        
           recvCount = localCount           

#else /* #ifndef INVERSE_GCELL_FILL */

        recvCount = 0
        if (gettingFrom .gt. 0) then

           bodyInfo % particles(FUL_PART_PROP,1:gettingFrom) = 0.
           bodyInfo % particles(FVL_PART_PROP,1:gettingFrom) = 0.
#if NDIM == MDIM
           bodyInfo % particles(FWL_PART_PROP,1:gettingFrom) = 0.
#endif           
           recvCount = gettingFrom
           localcount= recvCount
#endif  /* #ifndef INVERSE_GCELL_FILL */

           do p = 1, recvCount

#ifndef INVERSE_GCELL_FILL 
              blockID = int(particles2(BLK_PART_PROP,p))
              particleData =  particles2(1:NPART_PROPS,p)
#else
#endif

              !call Timers_start("Grid_updateSolidBodyForces") 
              call Grid_updateSolidBodyForces(b,p,blockID, particleData)
              !call Timers_stop("Grid_updateSolidBodyForces") 

#ifndef INVERSE_GCELL_FILL
              if (p .le. gettingFrom) &
              bodyInfo % particles(1:NPART_PROPS,p) = particleData 
#else
#endif

           enddo

        end if

     end if ! if (bodyInfo % myPE == bodyInfo % bodyMaster) then


  Particles_Allbods = Particles_Allbods + recvCount

  end do  ! End Loop over Bodies
! call Timers_stop("Force_particles")

! call Timers_start("Barrier_UpForces")
 lb = Particles_Allbods
 call MPI_Allreduce(lb,Particles_Allbods, CONSTANT_ONE, FLASH_INTEGER,MPI_MAX, MPI_COMM_WORLD,ierr)
! call Timers_stop("Barrier_UpForces")

 if(gr_meshMe .eq. MASTER_PE) write(*,*) 'MAX Particles per Proc=',Particles_Allbods !,FzTotL


#endif /* DO_MARKER_FORCING */

  ! Force inside Bodies:
!  call Timers_start("ib_forceInsideBody")
  call ib_forceInsideBody()
!  call Timers_stop("ib_forceInsideBody")


  !! Do the sum of forces to Velocities:
  !! ----------------------------------
  do lb=1,count

     blockID = listOfBlocks(lb)

     ! Get face data (velocities):
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 

     ! X velocities:
     facexData(VELC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1, &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &                       
     facexData(VELC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1, &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) + &  
     ib_dt*                                                                &
     facexData(FORC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1, &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))

     ! Y velocities:
     faceyData(VELC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &                       
     faceyData(VELC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) + &  
     ib_dt*                                                                &
     faceyData(FORC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))

#if NDIM == MDIM
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     ! Z velocities:
     facezData(VELC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1) = &                       
     facezData(VELC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1) + &  
     ib_dt*                                                                &
     facezData(FORC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1)
#endif

     ! Release face data (velocities):
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo

!  call Timers_stop("update_forces")
end subroutine gr_sbUpdateForces
