!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testI/sim_launchRays
!!
!!  NAME 
!!
!!   sim_launchRays
!!
!!  SYNOPSIS
!!
!!   sim_launchRays ()
!!
!!  DESCRIPTION
!!
!!   This routine launches 8 rays (2 rays in 2D) on a circle onto the xz-plane of a quadratic
!!   potential tube and records their (x,z) exit points and their exit power. The analysis of
!!   these results is done in a separate routine.
!!
!! ARGUMENTS
!!
!!***

subroutine sim_launchRays ()

  use Simulation_data,             ONLY : sim_globalComm,          &
                                          sim_globalMe,            &
                                          sim_nFocalRays,          &
                                          sim_nRaysMax,            &
                                          sim_rayPexit,            &
                                          sim_rayXexit,            &
                                          sim_rayZexit,            &
                                          sim_XZtypeRays,          &
                                          sim_yfocal

  use ed_extractRaysData,          ONLY : ed_extractRayData

  use Driver_interface,            ONLY : Driver_abortFlash

  use Grid_interface,              ONLY : Grid_getListOfBlocks

  use EnergyDeposition_interface,  ONLY : EnergyDeposition

  implicit none

#include "EnergyDeposition.h"
#include "Flash.h"
#include "constants.h"

  include "Flash_mpi.h"

  logical  :: success

  integer  :: blockCount
  integer  :: error
  integer  :: messageTag
  integer  :: nTimeSteps
  integer  :: pass
  integer  :: ray
  integer  :: request
  integer  :: timestep

  real     :: dt
  real     :: Pexit
  real     :: time
  real     :: Xexit, Yexit, Zexit

  integer  :: blockList (1:MAXBLOCKS)
  integer  :: status    (1:MPI_STATUS_SIZE)

  real     :: sendarray (1:3)
  real     :: recvarray (1:3)
!
!
!     ...Get the list of blocks on current processor.
!
!
  call Grid_getListOfBlocks (LEAF,   blockList, blockCount)
!
!
!     ...Call the energy deposition routine. Each time step activates one of the beams, each
!        of which launches 1 ray at a specific location:
!
!        For 2D cases: launch a total of 2 rays located on 3cm away from
!                      the x-axis midpoint 5cm. Hence the launching (x)
!                      locations are (in time order):
!
!                                    1) ---> (8)
!                                    2) ---> (2)
!
!        For 3D cases: launch a total of 8 rays located on the circle around
!                      the tube's xz-plane midpoint (5cm,5cm) and radius of 3cm.
!                      Detailed (x,z) locations of launching are (in time order):
!
!                              1) ---> (8,5)
!                              2) ---> (2,5)
!                              3) ---> (5,8)
!                              4) ---> (5,2)
!                              5) ---> (5 + 3/sqrt(2) , 5 + 3/sqrt(2))
!                              6) ---> (5 + 3/sqrt(2) , 5 - 3/sqrt(2))
!                              7) ---> (5 - 3/sqrt(2) , 5 + 3/sqrt(2))
!                              8) ---> (5 - 3/sqrt(2) , 5 - 3/sqrt(2))
!
!        During each time step:
!
!           1) At the beginning, the ray array is inititalized to zero on all processors,
!              and only on 1 processor the first ray array entry is active, containing
!              the initital data of the ray being launched.
!
!           2) At the end, some of the processors will have ray data in the first array
!              entry corresponding to the block exit points when the ray is travelling
!              throught the domain. We are interested only in the ray array entry which
!              corresponds to the true xz-plane domain exit point at y = yfocal. Only one
!              processor can have this specific xz-plane domain exit point.
!
!
  if (NDIM == 2) then
      nTimeSteps = 2
      sim_XZtypeRays = .false.
  end if

  if (NDIM == 3) then
      nTimeSteps = 8
      sim_XZtypeRays = .true.
  end if
  
  if (sim_globalMe == MASTER_PE) then
      sim_nFocalRays = 0
  end if

  dt   = 1.e-8      ! this ensures correct initial ray power of 1 erg
  pass = 1          ! mimics unsplit energy deposition
  time = 0.0        ! overall simulation time

  do timestep = 1,nTimeSteps

     call EnergyDeposition (blockCount, blockList, dt, time, pass)

     call ed_extractRayData (rayID      = 1,        &
                             entryField = RAY_POSX, &
                             dataValue  = Xexit,    &
                             dataFound  = success   )

     if (.not.success) then
          call Driver_abortFlash ('[sim_launchRays] ERROR: Could not extract rays X position')
     end if

     call ed_extractRayData (rayID      = 1,        &
                             entryField = RAY_POSY, &
                             dataValue  = Yexit,    &
                             dataFound  = success   )

     if (.not.success) then
          call Driver_abortFlash ('[sim_launchRays] ERROR: Could not extract rays Y position')
     end if

     call ed_extractRayData (rayID      = 1,        &
                             entryField = RAY_POSZ, &
                             dataValue  = Zexit,    &
                             dataFound  = success   )

     if (.not.success) then
          call Driver_abortFlash ('[sim_launchRays] ERROR: Could not extract rays Z position')
     end if

     call ed_extractRayData (rayID      = 1,        &
                             entryField = RAY_POWR, &
                             dataValue  = Pexit,    &
                             dataFound  = success   )

     if (.not.success) then
          call Driver_abortFlash ('[sim_launchRays] ERROR: Could not extract rays power')
     end if

     if (Yexit == sim_yfocal) then

         sendarray (1) = Xexit
         sendarray (2) = Zexit
         sendarray (3) = Pexit

         messageTag = 1

         call MPI_Isend  (sendarray,      &
                          3,              &
                          FLASH_REAL,     &
                          MASTER_PE,      &
                          messageTag,     &
                          sim_globalComm, &
                          request,        &
                          error           )
     end if

     if (sim_globalMe == MASTER_PE) then

         call MPI_Recv (recvarray,       &
                        3,               &
                        FLASH_REAL,      &
                        MPI_ANY_SOURCE,  &
                        MPI_ANY_TAG,     &
                        sim_globalComm,  &
                        status,          &
                        error            )

         sim_nFocalRays = sim_nFocalRays + 1
         
         if (sim_nFocalRays > sim_nRaysMax) then
             call Driver_abortFlash ('[sim_launchRays] ERROR: Cannot store all focal rays')
         end if

         sim_rayXexit (sim_nFocalRays) = recvarray (1)
         sim_rayZexit (sim_nFocalRays) = recvarray (2)
         sim_rayPexit (sim_nFocalRays) = recvarray (3)

     end if

     time = time + dt

     call MPI_Barrier (sim_globalComm , error)

  end do
!
!
!     ...Check, if all focal rays data is on the master processor. At this stage, only
!        the master processor has the information of all rays launched and in proper order.
!
!
  if (sim_globalMe == MASTER_PE) then

      if (sim_nFocalRays /= nTimeSteps) then
          call Driver_abortFlash ('[sim_launchRays] ERROR: Not all rays on master processor')
      end if

  end if
!
!
!     ...Ready!
!
!
  return
end subroutine sim_launchRays
