!!****if* source/Grid/GridSolvers/Pfft/RegularGridSolver/gr_pfftSpecifyTransform
!!
!! NAME
!!
!! gr_pfftSpecifyTransform
!!
!! SYNOPSIS
!!  
!! gr_pfftSpecifyTransform(integer(OUT) :: transformType(MDIM))
!! 
!! DESCRIPTION
!!
!! Initialises the transform type for the particular solver.
!!
!! ARGUMENTS
!!  
!! transformType - Array containing the type of transform required 
!!                 along each axis (e.g. PFFT_REAL2C or PFFT_COMPLEX)
!!
!! NOTES
!!
!! Called by gr_pfftInit.
!!
!!***

subroutine gr_pfftSpecifyTransform (transformType, baseDatType, bcTypes)

  use Grid_interface,   ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
                               GRID_PDE_BND_DIRICHLET
  use Grid_data,        ONLY : gr_solverType, gr_iSolveType, gr_jSolveType, gr_kSolveType, &
                               gr_iStr, gr_jStr, gr_kStr
  use gr_pfftData,      ONLY : pfft_solver
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Pfft.h"

  integer, dimension(MDIM), intent(OUT) :: transformType
  integer, dimension(0:MDIM), intent(OUT), OPTIONAL :: baseDatType
  integer, dimension(2*MDIM), intent(IN),  OPTIONAL :: bcTypes
  integer :: i

  if (.true.) then
  
    ! specify poisson solver combination
    select case (gr_solverType)

    ! automatic selection based on stretching
    case (AUTO)
      pfft_solver = TRIG_TRIG 
      if (gr_iStr) pfft_solver += 1
      if (gr_jStr) pfft_solver += 1
      if (NDIM .eq. 3) then
        pfft_solver += (TRIG_TRIG_TRIG - TRIG_TRIG)
        if (gr_kStr) pfft_solver += 1
      endif

    ! user specified selection
    case (USER_DEFINED)
      pfft_solver = TRIG_TRIG
      if (gr_iSolveType .eq. MATRIX) pfft_solver += 1
      if (gr_jSolveType .eq. MATRIX) pfft_solver += 1
      if (NDIM .eq. 3) then
        pfft_solver += (TRIG_TRIG_TRIG - TRIG_TRIG)
        if (gr_kSolveType .eq. MATRIX) pfft_solver += 1
      endif
    
    ! fully trignometric poisson solver 
    case (TRIG)
      pfft_solver = TRIG_TRIG
      if (NDIM .eq. 3) pfft_solver += (TRIG_TRIG_TRIG - TRIG_TRIG)

    ! either pentadiagonal or septadiagonal poisson solver
    case (MATRIX)
      pfft_solver = DRCT_DRCT
      if (NDIM .eq. 3) pfft_solver = DRCT_DRCT_DRCT

    ! unsupported method of specifying poisson solver
    case default
      pfft_solver = 0
      call Driver_abortFlash("This Poisson Solver requires poisson_solver == &
                             &auto, user, trig, or direct!")  
    end select

    ! verify supported poisson solver combination
    if (.not. (                                                                    &
       ((pfft_solver .ge. TRIG_TRIG)      .and. (pfft_solver .le. DRCT_DRCT)) .or. & 
       ((pfft_solver .ge. TRIG_TRIG_TRIG) .and. (pfft_solver .le. DRCT_DRCT_DRCT)) &
       )) then
      call Driver_abortFlash("Specification of Poisson Solver options not supported!")  
    endif
 
 endif
  
  if (.true.) then

    !! In this implementation we are only working with matching boundary conditions
    if (present(bcTypes)) then
      if (bcTypes(1) /= bcTypes(2)) then
        if ((bcTypes(1)==GRID_PDE_BND_NEUMANN.OR.bcTypes(1)==GRID_PDE_BND_DIRICHLET) .NEQV. &
            (bcTypes(2)==GRID_PDE_BND_NEUMANN.OR.bcTypes(2)==GRID_PDE_BND_DIRICHLET)) &
          call Driver_abortFlash("This Poisson solver requires the same type of boundaries left and right!")
      end if
      if (bcTypes(1) /= GRID_PDE_BND_PERIODIC .AND. bcTypes(1) /= GRID_PDE_BND_NEUMANN &
                                              .AND. bcTypes(1) /= GRID_PDE_BND_DIRICHLET) &
          call Driver_abortFlash(&
          "This Poisson solver requires periodic or homogeneous Dirichlet or Neumann boundaries in the X direction!")
#if NDIM > 1
      if (bcTypes(3) /= bcTypes(4)) then
        if ((bcTypes(3)==GRID_PDE_BND_NEUMANN.OR.bcTypes(3)==GRID_PDE_BND_DIRICHLET) .NEQV. &
            (bcTypes(4)==GRID_PDE_BND_NEUMANN.OR.bcTypes(4)==GRID_PDE_BND_DIRICHLET)) &
          call Driver_abortFlash("This Poisson solver requires the same type of boundaries up and down!")
      end if
      if (bcTypes(3) /= GRID_PDE_BND_PERIODIC .AND. bcTypes(3) /= GRID_PDE_BND_NEUMANN &
                                              .AND. bcTypes(3) /= GRID_PDE_BND_DIRICHLET) &
          call Driver_abortFlash(&
          "This Poisson solver requires periodic or homogeneous Dirichlet or Neumann boundaries in the Y direction!")
#endif
#if NDIM > 2
      if (bcTypes(5) /= bcTypes(6)) then
        if ((bcTypes(5)==GRID_PDE_BND_NEUMANN.OR.bcTypes(5)==GRID_PDE_BND_DIRICHLET) .NEQV. &
            (bcTypes(6)==GRID_PDE_BND_NEUMANN.OR.bcTypes(6)==GRID_PDE_BND_DIRICHLET)) &
          call Driver_abortFlash("This Poisson solver requires the same type of boundaries front and back!")
      end if
      if (bcTypes(5) /= GRID_PDE_BND_PERIODIC .AND. bcTypes(5) /= GRID_PDE_BND_NEUMANN &
                                            .AND. bcTypes(5) /= GRID_PDE_BND_DIRICHLET) &
          call Driver_abortFlash(&
          "This Poisson solver requires periodic or homogeneous Dirichlet or Neumann boundaries in the Z direction!")
#endif
      do i=1, NDIM
        select case (bcTypes(2*i - 1))
        case(GRID_PDE_BND_PERIODIC)
          transformType(i) = PFFT_REAL
        case(GRID_PDE_BND_NEUMANN)
          if (bcTypes(2*i)==GRID_PDE_BND_DIRICHLET) then
            transformType(i) = PFFT_COS_IV
          else
            transformType(i) = PFFT_COS_CC
          end if
        case(GRID_PDE_BND_DIRICHLET)
          if (bcTypes(2*i)==GRID_PDE_BND_NEUMANN) then
            transformType(i) = PFFT_SIN_IV
          else
            transformType(i) = PFFT_SIN_CC
          end if
        case default
          call Driver_abortFlash("This Poisson solver requires periodic or &
                                 &homogeneous Dirichlet or Neumann boundaries!")
        end select
      end do
    else
      transformType(IAXIS) = PFFT_REAL
      transformType(JAXIS:KAXIS) = PFFT_REAL
    end if
    if (present(baseDatType)) then
      baseDatType(0:MDIM) = PFFT_PCLDATA_REAL
    endif
  end if

end subroutine gr_pfftSpecifyTransform
