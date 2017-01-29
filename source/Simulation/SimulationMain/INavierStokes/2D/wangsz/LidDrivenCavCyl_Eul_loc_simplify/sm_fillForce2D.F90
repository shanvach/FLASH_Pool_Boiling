!!****if* source/Grid/GridMain/UG/Grid_fillGuardCells
!!
!! NAME
!!  Grid_fillGuardCells
!!
!! SYNOPSIS
!!
!!  call Grid_fillGuardCells(integer(IN) :: gridDataStruct,
!!                           integer(IN) :: idir,
!!                  optional,integer(IN) :: minLayers,
!!                  optional,integer(IN) :: eosMode,
!!                  optional,logical(IN) :: doEos,
!!                  optional,integer(IN) :: maskSize,
!!                  optional,logical(IN) :: mask(maskSize),
!!                  optional,logical(IN) :: makeMaskConsistent,
!!                  optional,integer(IN) :: selectBlockType,
!!                  optional,logical(IN) :: unitReadsMeshDataOnly)
!!
!!
!! DESCRIPTION 
!!  
!!  The argument "gridDataStruct" can take on one of many valid 
!!  values to determine a specific grid data structure on which to apply
!!  the guardcell fill operation. The currently available options are listed with
!!  the arguments. Most users will use CENTER as the option, 
!!  since applications typically use the cell centered grid data, and they want
!!  guardcells to be filled for all the variables.
!!  More specialized applications, such as the unsplit methods, may want to use
!!  other options. 
!!  The user can also choose to fill guard cells either in a single direction,
!!  or all of them. 
!!
!!
!! ARGUMENTS 
!!  
!!
!!  gridDataStruct - integer constant, defined in "constants.h", 
!!                   indicating which grid data structure 
!!                   variable's guardcells to fill.
!!                   UG has 4 data structures for grid variables that
!!                   can have their guardcells filled. 
!!
!!                   unk                all cell centered variables in the grid
!!                   facex,facey,facez  all face centered variables along i,j,k 
!!                                      direction respectively
!!                   
!!                   valid values of gridDataStruct are  
!!                   CENTER             unk only
!!                   WORK               has no meaning in UG
!!                   The face variables are not yet implemented in UG
!!                   FACES              facex,facey, and facez
!!                   FACEX              facex
!!                   FACEY              facey
!!                   FACEZ              facez
!!                   CENTER_FACES     unk,facex,facey,facez
!!
!!  idir - direction of guardcell fill.  User can specify ALLDIR for all (x,y,z)
!!         directions, or if for example the algorithm only does one directional
!!         sweep at a time then time can be saved by filling only the guardcell
!!         direction that is needed.  A user would pass in the constants defined
!!         in constants.h IAXIS, JAXIS or KAXIS to fill guardcells in only one 
!!         direction.        
!!         All layers of guardcells in the given direction(s) are filled.
!!         In the current UG implementation, idir is ignored and a full
!!         guardcell fill in all directions is always performed.
!!
!!         THE REMAINING ARGUMENTS HAVE NO MEANING IN UG
!!
!!  minLayers - number of guardcell layers requested for all directions.
!!
!!   eosMode  - The mode in which eos is to be applied
!!   doEos    - the UG implementation does not act upon this argument
!!   maskSize - the size of the mask array. 
!! 
!!  mask -  It is a one-dimensional logical array 
!!          with indices corresponding to variables in the grid data
!!          structures. If a variable should have its guardcells filled,
!!          the corresponding element in "mask" is true, otherwise it is
!!          false.
!!          The mask is always ignored if the runtime parameter
!!          enableMaskedGCFill is set .FALSE.
!!  
!! makeMaskConsistent - If true when mask is applied, it is made sure that for
!!          all the selected variables in the mask, the ones they are dependent
!!          on are true too. It is also determined whether there is a need to 
!!          apply Eos if doEos argument is true.
!!
!! selectBlockType - IGNORED
!!
!! unitReadsMeshDataOnly - specifies that the unit calling Grid_fillGuardCells
!!                         does not update any internal grid data.  This
!!                         allows us to skip the next guard cell fill because
!!                         the guard cells already contain up to date data.
!!
!! EXAMPLE
!!
!!   #include "Flash.h"
!!   #include "constants.h"
!!
!!      call Grid_fillGuardCells( CENTER, IAXIS)
!!
!!     This call will fill all guardcells for all cell-centered 
!!     variables in the x direction.
!!     
!! EXAMPLE 2
!!
!!   #include "Flash.h"
!!   #include "constants.h"
!!
!!      call Grid_fillGuardCells( CENTER_FACES, ALLDIR)
!!     
!!     This call fills guardcells along all directions in both
!!     cell centered and face centered data structures.
!!
!! NOTES
!!
!!   The masking functionality is not yet included in UG
!!  
!!***


#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine sm_fillForce2D(fX, fY, tmpX, tmpY)

  use Grid_data, ONLY : gr_axisComm, gr_exch, gr_exchLoc, gr_gridDataStruct, &
       gr_justExchangedGC,gr_domainBC, &
       gr_offset,gr_allPeriodic,gr_bndOrder, gr_meshMe
  use Grid_interface, ONLY : Grid_getBlkIndexLimits
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_bcInterface, ONLY : gr_bcApplyToAllBlks

  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none
#include "constants.h"
#include "Flash.h"

  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC), intent(in) :: fX
  real, dimension(GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC), intent(in) :: fY
  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC), intent(out) :: tmpX
  real, dimension(GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC), intent(out) :: tmpY

  ! integer that hold the boundary condition type (ie, PERIODIC)
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: blockID = 1
  
  !These are arrays to hold the sending start point indicies and receiving
  !start point indicies for the MPI data types in the shift data subroutine
  !They are used to specify the starting incidies in the unk array
  !example 
  !sendRight(1,1) in x direction holds the 1st unk dimension which is the variable ID
  !sendRight(1,2) in x direction holds the 2nd unk dimension which is the i starting index
  !sendRight(1,3) in x direction holds the 3rd unk dimension which is the j starting index
  !sendRight(1,4) in x direction holds the 4th unk dimension which is the k starting index
  integer, dimension(MDIM, MDIM+1) :: sendRight, sendLeft, recvRight, recvLeft


  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,FACEX)
  recvLeft(:,:) = 1
  sendLeft(:,:) = 1
  recvRight(:,:) = 1
  sendRight(:,:) = 1

  sendRight(2,3) = blkLimits(HIGH,JAXIS) + 1
  sendLeft(2,3)  = blkLimits(LOW,JAXIS) - 2

  recvRight(2,3)  = blkLimits(LOW,JAXIS)
  recvLeft(2,3)   = blkLimits(HIGH,JAXIS) - 1

  call gr_shiftData_loc(gr_axisComm(JAXIS), gr_exchLoc(FACEX_DATATYPE,JAXIS), & 
                  &     sendRight(JAXIS,:), sendLeft(JAXIS,:),  & 
                  &     recvRight(JAXIS,:), recvLeft(JAXIS,:),  & 
                  &     GRID_IHI_GC+1, GRID_JHI_GC, GRID_KHI_GC, fX, tmpX)

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,FACEY)
  recvLeft(:,:) = 1
  sendLeft(:,:) = 1
  recvRight(:,:) = 1
  sendRight(:,:) = 1

  sendRight(2,3) = blkLimits(HIGH,JAXIS) 
  sendLeft(2,3)  = blkLimits(LOW,JAXIS) - 1

  recvRight(2,3)  = blkLimits(LOW,JAXIS)
  recvLeft(2,3)   = blkLimits(HIGH,JAXIS) - 1

  call gr_shiftData_loc(gr_axisComm(JAXIS), gr_exchLoc(FACEY_DATATYPE,JAXIS), & 
                  &     sendRight(JAXIS,:), sendLeft(JAXIS,:),  & 
                  &     recvRight(JAXIS,:), recvLeft(JAXIS,:),  & 
                  &     GRID_IHI_GC, GRID_JHI_GC+1, GRID_KHI_GC, fY,tmpY)


  if(.not.gr_allPeriodic) then
    call gr_bcApplyToAllBlks_loc2D(tmpX, tmpY)
  end if

  return
end subroutine sm_fillForce2D


