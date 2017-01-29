!!***if* source/Simulation/SimulationMain/PhoenixInput/Simulation_finalize
!!
!! NAME
!!  Simulation_finalize
!!
!! SYNOPSIS
!!
!!  Simulation_finalize()
!!
!! DESCRIPTION
!!
!!  This dummy function cleans up the Simulation unit, deallocates memory, etc.
!!  However, as nothing needs to be done, only this stub is included.
!!
!!  In this case it is freeing up the property names from the nucleosynthetic 
!!  file.
!!
!! ARGUMENTS
!!
!!
!!
!!***

subroutine Simulation_finalize()

  use Logfile_interface, ONLY : Logfile_stamp
  use Simulation_data, ONLY : sim_propNames

  implicit none

  integer :: istat

  deallocate(sim_propNames,stat=istat)
  if (istat .NE. 0) then
     call Logfile_stamp('WARNING - failed to deallocate sim_propNames','[Simulation_finalize]')
  end if
  

  return

end subroutine Simulation_finalize
