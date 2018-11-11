!!****f* source/physics/IncompNS/IncompNS_finalize
!!
!! NAME
!!
!!  IncompNS_finalize
!!
!!
!! SYNOPSIS
!!
!!  IncompNS_finalize()
!!  
!!
!! DESCRIPTION
!! 
!!  Finalize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_finalizeFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***

subroutine IncompNS_finalize()

  use Multiphase_data, only: mph_timeStampAll, mph_nucSiteTemp, mph_isAttachedOld, mph_isAttachedAll

  implicit none

  deallocate(mph_timeStampAll)
  deallocate(mph_nucSiteTemp)
  deallocate(mph_isAttachedOld)
  deallocate(mph_isAttachedAll)

end subroutine IncompNS_finalize

