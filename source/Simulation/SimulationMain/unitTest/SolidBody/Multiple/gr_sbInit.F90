!!****if* source/Simulation/SimulationMain/unitTest/SolidBody/Multiple/gr_sbInit
!!
!! NAME
!!  gr_sbInit
!!
!!
!! SYNOPSIS
!!  gr_sbInit()
!!
!!
!! DESCRIPTION
!!
!! Called from Grid_init. Randomize the initial particle positions of each
!! solid body.  Allocate and populate the data structure that holds all
!! information about each solid body. 
!!
!! ARGUMENTS      
!!
!! PARAMETERS
!!
!!***

!Called from Grid_init
!Allocate and populate the data structure that holds all information
!about the single solid body.

#include "constants.h"
#include "Flash.h"

Subroutine gr_sbInit()
  use Grid_data, ONLY : gr_meshComm
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, &
       gr_sbPtNumX, gr_sbPtNumY, gr_sbPtNumZ, gr_sbDebug
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none
  include "Flash_mpi.h"

  integer :: ierr, id, b
  real :: x, y, z, x_min, x_max, y_min, y_max, z_min, z_max, x_min_final, &
       x_max_final, y_min_final, y_max_final, z_min_final, z_max_final

  call RuntimeParameters_get("sb_NumBodies", gr_sbNumBodies)
  allocate(gr_sbBodyInfo(gr_sbNumBodies))
  call random_seed()
  do b = 1, gr_sbNumBodies
     gr_sbBodyInfo(b) % boundBox(:,:) = 0.0
     call MPI_Comm_rank(gr_meshComm, id, ierr)
     if (id == 0) then ! randomize particle positions
        call random_number(x)
!        print *, "body number", b, "random number", x
        x_min = 0.1+(0.9-0.1)*x ! random number between 0.1 and 0.9
        x_min_final = (int(10*x_min))/10.0 
        if (x_min_final == 0.9) then
           x_max_final = 1.0
        else
          x_max = (x_min+0.1) + (1-(x_min+0.1))*x
           x_max_final = (int(x_max*10))/10.0
        end if
        call random_number(y)
        y_min = 0.1+(0.9-0.1)*y
        y_min_final = (int(10*y_min))/10.0
        if (y_min_final == 0.9) then
           y_max_final = 1.0
        else
           y_max = (y_min+0.1) + (1-(y_min+0.1))*y
           y_max_final = (int(y_max*10))/10.0
        endif
        call random_number(z)
        z_min = 0.1+(0.9-0.1)*z
        z_min_final = (int(10*z_min))/10.0
        if (z_min_final == 0.9) then
           z_max_final = 1.0
        else
           z_max = (z_min+0.1) + (1-(z_min+0.1))*z
           z_max_final = (int(z_max*10))/10.0
        endif
     endif
     call MPI_BCAST(x_min_final, 1, FLASH_REAL, 0, gr_meshComm, ierr)
     call MPI_BCAST(x_max_final, 1, FLASH_REAL, 0, gr_meshComm, ierr)
     call MPI_BCAST(y_min_final, 1, FLASH_REAL, 0, gr_meshComm, ierr)
     call MPI_BCAST(y_max_final, 1, FLASH_REAL, 0, gr_meshComm, ierr)
     call MPI_BCAST(z_min_final, 1, FLASH_REAL, 0, gr_meshComm, ierr)
     call MPI_BCAST(z_max_final, 1, FLASH_REAL, 0, gr_meshComm, ierr)

     gr_sbBodyInfo(b) % boundBox(LOW,IAXIS) = x_min_final
     gr_sbBodyInfo(b) % boundBox(HIGH,IAXIS) = x_max_final
     call RuntimeParameters_get("sb_ptNumX", gr_sbPtNumX)
   
     if (NDIM >= 2) then
        gr_sbBodyInfo(b) % boundBox(LOW,JAXIS) = y_min_final
        gr_sbBodyInfo(b) % boundBox(HIGH,JAXIS) = y_max_final
        call RuntimeParameters_get("sb_ptNumY", gr_sbPtNumY)
     else
        gr_sbPtNumY = 1
     end if

     if (NDIM == 3) then
        gr_sbBodyInfo(b) % boundBox(LOW,KAXIS) = z_min_final
        gr_sbBodyInfo(b) % boundBox(HIGH,KAXIS) = z_max_final
        call RuntimeParameters_get("sb_ptNumZ", gr_sbPtNumZ)
     else
        gr_sbPtNumZ = 1
     end if

     call RuntimeParameters_get("sb_debug", gr_sbDebug)
  enddo
End Subroutine gr_sbInit
