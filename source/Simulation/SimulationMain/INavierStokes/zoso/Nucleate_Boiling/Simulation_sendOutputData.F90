!!****f* source/Simulation/Simulation_sendOutputData
!!
!! NAME
!!  Simulation_sendOutputData
!!
!! SYNOPSIS
!! 
!!  Simulation_sendOutputData()
!!  
!! DESCRIPTION 
!!
!! This routine sends the scalar variables owned by the Simulation unit
!! to the IO unit, to be written to a checkpoint file.
!!
!!
!!***

subroutine Simulation_sendOutputData()

    use Heat_AD_data, only: ht_fmic,ht_qmic
    use IO_interface, ONLY :  IO_setScalar

    implicit none

    call IO_setScalar("qmic", ht_qmic)
    call IO_setScalar("fmic", ht_fmic)

end subroutine Simulation_sendOutputData
