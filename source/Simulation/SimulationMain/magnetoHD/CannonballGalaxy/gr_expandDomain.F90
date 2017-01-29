!!****if* source/Simulation/SimulationMain/magnetoHD/CannonballGalaxy/gr_expandDomain
!!
!!  NAME
!!     gr_expandDomain
!!
!!  SYNOPSIS
!!     call gr_expandDomain(logical,(OUT) :: particlesInitialized)
!!
!!  DESCRIPTION
!!
!!    The grid is initialized in gr_createDomain, with specified
!!    number of blocks. Typically single block. This routine
!!    refines appropriate portions of the initialized physical 
!!    domain according to the given criterion, and applies initial 
!!    conditions to the AMR domain. 
!!
!!  ARGUMENTS
!!    particlesInitialized : is true if this routine initialized particles positions
!!
!!***

!!REORDER(4):solnData

#define DEBUG_PARTICLES

subroutine gr_expandDomain (particlesInitialized)


#include "Flash.h"
  use Grid_data, ONLY : gr_domainBC,gr_eosModeInit,gr_refineOnParticleCount,&
       gr_refineOnPdens,gr_maxParticlesPerBlk,gr_minParticlesPerBlk, gr_meshMe,&
       gr_meshNumProcs, gr_lrefineMinInit
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_stampVarMask
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getLocalNumBlks, Grid_getListOfBlocks, Grid_fillGuardCells, &
       Grid_markRefineDerefine, Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getDeltas, Grid_getBlkPhysicalSize, Grid_getBlkCenterCoords
  use Grid_data, ONLY : gr_ihiGc,gr_jhiGc,gr_khiGc,gr_blkList
  use tree, ONLY : lrefine, lrefine_min, lrefine_max, grid_changed
  use paramesh_interfaces, ONLY : amr_refine_derefine, amr_restrict
  use Eos_interface, ONLY : Eos_wrapped
  use Simulation_data, ONLY : sim_plasmaBeta, sim_xCtr, sim_yCtr, sim_zCtr

#ifndef FLASH_GRID_PARAMESH2
  use physicaldata, ONLY: no_permanent_guardcells, mpi_pattern_id
#endif
  use Simulation_interface, ONLY : Simulation_initBlock, Simulation_customizeProlong
  use Particles_interface, ONLY : Particles_accumCount, &
    Particles_initPositions, &
    Particles_updateRefinement
  use Driver_interface, ONLY : Driver_abortFlash
  use RadTrans_interface, ONLY: RadTrans_sumEnergy
  implicit none

#include "constants.h"

  include 'Flash_mpi.h'

  real, pointer:: solnData(:,:,:,:), facexData(:,:,:,:), faceyData(:,:,:,:), &
       facezData(:,:,:,:)
  logical, intent(out) :: particlesInitialized
  integer :: lnblocks, lrefineMinSave,j,k,l,blockID

  !!          Local variables and functions

  integer :: ntimes, i, ierr

  integer, dimension(2,MDIM) :: blkLimits
  integer, dimension(2,MDIM) :: blkLimitsGC
  integer :: count, cur_treedepth, grid_changed_anytime
  logical :: restart = .false.
  logical :: particlesPosnsDone, updateRefine
  integer :: level = FINEST  !not yet implemented, 1 is a dummy value
  integer ,dimension(MAXBLOCKS) :: blkList
  character(len=32), dimension(2,2) :: block_buff
  character(len=32)                 :: int_to_str
  integer :: gridDataStruct, whichBlocks
  real, dimension(MDIM) :: blockSize, del, blockCoord
  real :: totmagp, totpres, ltotmagp, ltotpres, rb

  !!============================================================================



  !!============================================================================

  !! The beginning timestep number, time, and timestep.
  !! If the initial redshift (zinitial) is physical (>= 0),
  !! we use it to initialize the time; otherwise we set the
  !! redshift to zero and get the initial time from tinitial.
  !! The latter case (no cosmology) is the default.

  ! initialize the step counter and the simulation time
  ! the timestep initialization is moved to after the initialization,
  ! so we can check whether it is > t_cfl

  particlesInitialized=.false.
  call Grid_getBlkIndexLimits(1, blkLimits, blkLimitsGC)

  call gr_initParameshArrays(restart,        &
       gr_domainBC(LOW,IAXIS),gr_domainBC(HIGH,IAXIS), &
       gr_domainBC(LOW,JAXIS),gr_domainBC(HIGH,JAXIS), &
       gr_domainBC(LOW,KAXIS),gr_domainBC(HIGH,KAXIS))

  ! The Paramesh call above may have already resulted in some block refining,
  ! so try get the current max level from the lrefine array. This is only used for
  ! diagnostic output. Note that this assumes that lrefine on the current 
  ! processor is representative of the grid as a whole.
  cur_treedepth = maxval(lrefine)

  gridDataStruct = CENTER
#if NFACE_VARS > 0
  gridDataStruct = CENTER_FACES
#endif

  lrefineMinSave = lrefine_min
  lrefine_min = min(gr_lrefineMinInit,lrefine_max)

  grid_changed_anytime = grid_changed ! save value established by previous Paramesh initialization

  updateRefine=.false.
  
  do ntimes = 1, lrefine_max+2
     if (ntimes .EQ. gr_lrefineMinInit) then
        lrefine_min = lrefineMinSave
     end if
     write (block_buff(1,1), '(a)') 'iteration'
     write (int_to_str, '(i7,a1)') ntimes, ','
     write (block_buff(1,2), '(a,1x,a)') trim(adjustl(int_to_str))
     
     write (block_buff(2,1), '(a)') 'create level'
     write (int_to_str, '(i7)') min(cur_treedepth+1,lrefine_max)
     write (block_buff(2,2), '(a)') trim(adjustl(int_to_str))

     call Logfile_stamp( block_buff, 2, 2, '[GRID gr_expandDomain]')

     call gr_updateData()
     call Grid_getLocalNumBlks(lnblocks)
     whichBlocks = LEAF
     ! Paramesh may have already refined the original root block(s) once by this point.
     ! (That can happen in particular when lrefine_min > 1.)
     ! So initialize all existing blocks in the first iteration here.  This makes sure
     ! that root blocks are not left with uninitialized contents. (Normally that
     ! situation would only last until Grid_fillGuardCells is called with LEAF blocks
     ! at refinement level 2 anyway, since PARENT blocks are then updated as a side
     ! effect by restriction.) - KW
     if (ntimes == 1) whichBlocks = ALL_BLKS
     call Grid_getListOfBlocks(whichBlocks, blkList,count)


#ifndef FLASH_GRID_PARAMESH2
     if (no_permanent_guardcells) then
        call gr_commSetup(gridDataStruct)
     end if
#endif

     do i = 1, count
        !  We need to zero data in case we reuse blocks from previous levels
        !  but don't initialize all data in Simulation_initBlock... in particular
        !  the total vs. internal energies can cause problems in the eos call that 
        !  follows.
        call Grid_getBlkPtr(blkList(i), solnData)
        solnData = 0.0
        call Grid_releaseBlkPtr(blkList(i), solnData)
        !      Now reinitialize the solution on the new grid so that it has
        !      the exact solution.
        call Simulation_initBlock (blkList(i))
     end do

#ifdef ERAD_VAR
     ! Sum radiation energy density over all meshes. This call is
     ! needed for mesh replication.
     call RadTrans_sumEnergy(ERAD_VAR, count, blkList)
#endif

     ! This is here for safety, in case the user did not take care to make things
     ! thermodynamically consistent in the initial state.- KW
     call Timers_start("eos")
     do i = 1, count
        call Eos_wrapped(gr_eosModeInit,blkLimits,blkList(i))
     end do

     call Timers_stop("eos")
!!     particlesPosnsDone=.false.
     if(gr_refineOnParticleCount ) then
        
        !!   This loop initializes the particle positions if
        !!   their count is one of the refinement criteria.
        !!   If the initialization routine intends to keep
        !!   the already initialized particles around, instead
        !!   of reinitializing them as the grid is refined,
        !!   it should return 
        !!   updateRefine true. If the whole file has been
        !!   read in then particlesPosnsDone should be true, otherwise false
        !!    
!!        particlesPosnsDone=particlesPosnsDone.and.updateRefine
        if(.not.updateRefine) then
           particlesPosnsDone=.false.
        end if
        call Particles_initPositions(particlesPosnsDone,updateRefine)
#ifdef DEBUG_PARTICLES
        if (gr_meshMe == MASTER_PE .OR. gr_meshNumProcs .LE. 4) then
           print*,'gr_expandDomain after Particles_initPositions on',gr_meshMe,':',particlesPosnsDone,updateRefine
        end if
#endif
     end if
     if (ntimes .le. lrefine_max+1) then
        ! Guard cell filling and Eos_wrapped are done in Grid_markRefineDerefine as needed.
        call Grid_markRefineDerefine()
        grid_changed_anytime = max(grid_changed, grid_changed_anytime)
        grid_changed = 0              ! will be 1 after amr_refine_derefine if the grid actually changed  
        call amr_refine_derefine()

        ! update the grid coordinates for the new mesh
        call Timers_start("updateData")
        call gr_updateData()
        call Timers_stop("updateData")

        call Simulation_customizeProlong(BEFORE)
        call amr_prolong (gr_meshMe, 1, NGUARD)
        call Simulation_customizeProlong(AFTER)

#ifndef FLASH_GRID_PARAMESH2
        if (grid_changed .NE. 0) mpi_pattern_id = -abs(mpi_pattern_id) !make it different from recognized values
#endif           
        if(gr_refineOnParticleCount.and.updateRefine) call Particles_updateRefinement(lnblocks)
        cur_treedepth = max(maxval(lrefine),min(cur_treedepth+1,lrefine_max))
        
     end if

  end do !ntimes

  grid_changed = max(grid_changed, grid_changed_anytime) !leave global flag true if grid changed in ANY iteration

  if(gr_refineOnParticleCount) then
     if(.not.particlesPosnsDone) call Driver_abortFlash(&
       "This distribution of particles will not fit on the grid,increase max_particles_per_blk, or decrease the particle count")
     particlesInitialized=.true.
  end if

  lrefine_min = lrefineMinSave

  call gr_ensureValidNeighborInfo(10)

#if 0
  call Grid_getListOfBlocks(LEAF,blkList,count)
  ltotmagp = 0.0
  ltotpres = 0.0
  do l = 1, count
     call Grid_getBlkPtr(blkList(l), solnData)
     call Grid_getBlkIndexLimits(blkList(l), blkLimits, blkLimitsGC)
     call Grid_getBlkPhysicalSize(blkList(l), blockSize)
     call Grid_getBlkCenterCoords(blkList(l), blockCoord)

     rb = sqrt((blockCoord(1)-sim_xCtr(1))**2 + &
          (blockCoord(2)-sim_yCtr(1))**2 + &
          (blockCoord(3)-sim_zCtr(1))**2)
          
     if (rb <= 0.75*sim_aC1) then

        ltotmagp = ltotmagp + sum(solnData(MAGP_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))) * &
             blockSize(1)*blockSize(2)*blockSize(3)
        ltotpres = ltotpres + sum(solnData(PRES_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &           
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))) * &
             blockSize(1)*blockSize(2)*blockSize(3)
        
     endif

     call Grid_releaseBlkPtr(blkList(l), solnData)
  enddo
  call MPI_Allreduce(ltotmagp, totmagp, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(ltotpres, totpres, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  print *, totmagp
  print *, totpres
  
  do l = 1, count
     blockID = blkList(l)
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
     call Grid_getBlkPtr(blockID,solnData)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
     call Grid_getDeltas(blockID,del)
     facexData(MAG_FACE_VAR,:,:,:) = facexData(MAG_FACE_VAR,:,:,:) / &
          sqrt(sim_plasmaBeta*totmagp/totpres)
     faceyData(MAG_FACE_VAR,:,:,:) = faceyData(MAG_FACE_VAR,:,:,:) / &
          sqrt(sim_plasmaBeta*totmagp/totpres)
     facezData(MAG_FACE_VAR,:,:,:) = facezData(MAG_FACE_VAR,:,:,:) / &
          sqrt(sim_plasmaBeta*totmagp/totpres)
     
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
               
              solnData(MAGX_VAR,i,j,k) = 0.5*(facexData(MAG_FACE_VAR,i+1,j,k) + &
                   facexData(MAG_FACE_VAR,i,j,k))
              solnData(MAGY_VAR,i,j,k) = 0.5*(faceyData(MAG_FACE_VAR,i,j+1,k) + &
                   faceyData(MAG_FACE_VAR,i,j,k))
              solnData(MAGZ_VAR,i,j,k) = 0.5*(facezData(MAG_FACE_VAR,i,j,k+1) + &
                   facezData(MAG_FACE_VAR,i,j,k))

           enddo
        enddo
     enddo
           
     solnData(MAGP_VAR,:,:,:) = 0.5*(solnData(MAGX_VAR,:,:,:)**2 + &
          solnData(MAGY_VAR,:,:,:)**2 + &
          solnData(MAGZ_VAR,:,:,:)**2)

     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

              solnData(DIVB_VAR,i,j,k) = (facexData(MAG_FACE_VAR,i+1,j,k) - &
                   facexData(MAG_FACE_VAR,i,j,k))/del(1) &
                   + (faceyData(MAG_FACE_VAR,i,j+1,k) - &
                   faceyData(MAG_FACE_VAR,i,j,k))/del(2) &
                   + (facezData(MAG_FACE_VAR,i,j,k+1) - &
                   facezData(MAG_FACE_VAR,i,j,k))/del(3)

           enddo

        enddo

     enddo

     call Grid_releaseBlkPtr(blockID,solnData)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
  enddo

#endif

  return
end subroutine gr_expandDomain
