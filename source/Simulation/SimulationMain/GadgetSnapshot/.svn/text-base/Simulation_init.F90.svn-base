!!****if* source/Simulation/SimulationMain/GadgetSnapshot/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  call Simulation_init()
!!
!!
!! DESCRIPTION
!!  Initializes all the parameters needed for a particular simulation.
!!  This version reads information and particle setup from a Gadget
!!  snapshot. For now this only works for N-body only (collisionless) 
!!  simulations. 
!!
!! ARGUMENTS
!!  
!!  none
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use Driver_interface, ONLY : Driver_getMype, Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_set
  use Cosmology_data, ONLY: csm_scaleFactor, csm_hubble, csm_omega, csm_lambda
  use Cosmology_interface, ONLY: Cosmology_redshiftToTime
  use Driver_data, ONLY: dr_initialSimTime, dr_simTime, dr_redshift, &
       dr_redshiftInitial

  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  

  integer :: success

  call Driver_getMype(MESH_COMM, sim_meshMe)

  ! Set translations from Gadget units

  sim_unitLength   = 3.085678e21
  sim_unitMass     = 1.989e43
  sim_unitVelocity = 1.0e5

  sim_unitTime     = sim_unitLength / sim_unitVelocity
  sim_unitDensity  = sim_unitMass * sim_unitLength**(-3.0)
  sim_unitPressure = sim_unitMass / sim_unitLength * &
       sim_unitTime**(-2.0)
  sim_unitEnergy   = sim_unitMass * sim_unitLength**2.0 * &
       sim_unitTime**(-2.0)

  ! Figure out what exactly we're reading in

  call RuntimeParameters_get('snapshotNumber', sim_snapshotNumber)
  call RuntimeParameters_get('numFiles', sim_numFiles)
  call RuntimeParameters_get('path', sim_path)
  call RuntimeParameters_get('basename', sim_basename)

  ! Read the Gadget header

  if (sim_meshMe == MASTER_PE) then
     call sim_readGadgetHeader(success, sim_snapshotNumber, sim_numFiles, &
          sim_path, sim_basename, sim_numParticles, sim_time, sim_redshift, &
          sim_hubble, sim_omegaMatter, sim_omegaLambda)

     if (success < 0) &
          call Driver_abortFlash("[Simulation_init]: Problem with reading Gadget header!")

     print *, "Number of particles = ", sim_numParticles       ! if MASTER_PE
     print *, "Initial simulation time = ", sim_time           ! if MASTER_PE
     print *, "Initial simulation redshift = ", sim_redshift   ! if MASTER_PE
     print *, "Hubble parameter = ", sim_hubble                ! if MASTER_PE
     print *, "Omega Matter = ", sim_omegaMatter               ! if MASTER_PE
     print *, "Omega Lambda = ", sim_omegaLambda               ! if MASTER_PE

  endif

  call MPI_Bcast(sim_numParticles, 1, FLASH_INT, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sim_time, 1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sim_redshift, 1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sim_hubble, 1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sim_omegaMatter, 1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sim_omegaLambda, 1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)

  ! If the snapshot gives us a box size, set it. If not, take it from the
  ! parameter database.

  if (sim_boxSize > 0.0) then
     sim_xMin = -0.5*sim_boxSize
     sim_xMax = 0.5*sim_boxSize
     sim_yMin = -0.5*sim_boxSize
     sim_yMax = 0.5*sim_boxSize
     sim_zMin = -0.5*sim_boxSize
     sim_zMax = 0.5*sim_boxSize
     call RuntimeParameters_set("xmin", sim_xMin)
     call RuntimeParameters_set("xmax", sim_xMax)
     call RuntimeParameters_set("ymin", sim_yMin)
     call RuntimeParameters_set("ymax", sim_yMax)
     call RuntimeParameters_set("zmin", sim_zMin)
     call RuntimeParameters_set("zmax", sim_zMax)
  else
     call RuntimeParameters_get("xmin", sim_xMin)
     call RuntimeParameters_get("xmax", sim_xMax)
     call RuntimeParameters_get("ymin", sim_yMin)
     call RuntimeParameters_get("ymax", sim_yMax)
     call RuntimeParameters_get("zmin", sim_zMin)
     call RuntimeParameters_get("zmax", sim_zMax)
  endif

  ! If we're using cosmology, then set the values of the parameters to whatever
  ! the snapshot had.

#ifdef COSMOLOGY
  call RuntimeParameters_get("useCosmology", sim_useCosmo)
#endif

  if (sim_useCosmo) then
     
     csm_hubble = sim_hubble * 100. * sim_unitVelocity / &
          (sim_unitLength * 1000.)
     csm_omega  = sim_omegaMatter
     csm_lambda = sim_omegaLambda
     
     call RuntimeParameters_set("zInitial",sim_redshift)
     
     csm_scaleFactor = 1./(1.+sim_redshift)
     dr_redshift = sim_redshift
     dr_redshiftInitial = sim_redshift

     call Cosmology_redshiftToTime(sim_redshift, sim_time)

  else
     
     sim_time = sim_time * sim_timeUnit
     
  endif

  ! Set the initial simulation time from the snapshot.

  dr_initialSimTime = sim_time
  dr_simTime        = sim_time

  return

end subroutine Simulation_init







