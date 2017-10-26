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

    use Heat_AD_data, only: ht_qmic,ht_fmic,ht_dxmin

    use Multiphase_data, only: mph_baseRadius, mph_baseCountAll, mph_isAttachedAll, mph_isAttachedOld, mph_timeStampAll

    use IO_interface, ONLY :  IO_setScalar

    implicit none

      call IO_setScalar("qmic", ht_qmic)
      call IO_setScalar("fmic", ht_fmic)
      call IO_setScalar("microdx",ht_dxmin)
      call IO_setScalar("baseradius",mph_baseRadius)
      call IO_setScalar("basecount",mph_baseCountAll)


end subroutine Simulation_sendOutputData

