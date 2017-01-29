!!****if* source/physics/ImBound/ImBoundMain/LagForce/parallel/ib_forceInsideBody
!!
!!
!! NAME
!!
!!  ib_forceInsideBody
!!
!!
!! SYNOPSIS
!!
!!  ib_forceInsideBody()
!!
!!
!! DESCRIPTION
!!
!! Routine that forces eulerian velocities inside bodies. Stub.
!! For the moment only analytical shapes are used.
!!
!!***

subroutine ib_forceInsideBody()

  use SolidMechanics_data, only : sm_MeshMe,sm_NumProcs, sm_NumBodies, sm_BodyInfo, sm_meshComm

  use Grid_interface, ONLY : Grid_getListOfBlocks,   &
                             Grid_getDeltas,         &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_getBlkBoundBox, Grid_getBlkCenterCoords

  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, only : Driver_abortFlash

  use ImBound_data, only : ib_dt, ib_BlockMarker_flag

  implicit none


  return

end subroutine ib_forceInsideBody
