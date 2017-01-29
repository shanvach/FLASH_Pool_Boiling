!!****if* source/Grid/GridMain/Grid_getMinCellSize
!!
!! NAME
!!  Grid_getMinCellSize
!!
!! SYNOPSIS
!!
!!  Grid_getMinCellSize(real (OUT)  :: minCellSize)
!!               
!!  
!! DESCRIPTION 
!!
!!  Returns the smallest possible cell size in a simulation in any dimension
!!  that does not represent an angle in curvilinear coordinates.
!!
!!
!! ARGUMENTS
!!
!!  minCellSize - returned value
!!
!!***

#include "Flash.h"

subroutine Grid_getMinCellSize(minCellSize)

  use Grid_data, ONLY : gr_minCellSize, gr_meshMe
  use Grid_interface,    ONLY : Grid_getLocalNumBlks

#ifdef FLASH_GRID_PARAMESH
  use tree, only : nodetype,lrefine,lrefine_max
#endif
  implicit none
#include "constants.h"
#include "Flash_mpi.h"

  real, intent(OUT) :: minCellSize


  integer            :: lb, i, mesh_lrefmax, lmax, ierr, lnblocks2

  mesh_lrefmax = 1  

#ifdef FLASH_GRID_PARAMESH
  call Grid_getLocalNumBlks(lnblocks2)
  ! Determine the largest level of refinement in the mesh.
  lmax = -1
  do i = 1, lnblocks2
  if (nodetype(i) == 1) then
    lmax = max( lrefine(i), lmax )
  endif
  enddo
  
  call mpi_allreduce ( lmax, mesh_lrefmax, 1, MPI_INTEGER, MPI_MAX, &
                       MPI_COMM_WORLD, ierr )

  minCellSize = real(lrefine_max)/real(mesh_lrefmax)*gr_minCellSize
  
  if (gr_meshMe .eq. MASTER_PE) write(*,*) "Min size=",lrefine_max,mesh_lrefmax,minCellSize

#else
  minCellSize = gr_minCellSize
#endif


end subroutine Grid_getMinCellSize
