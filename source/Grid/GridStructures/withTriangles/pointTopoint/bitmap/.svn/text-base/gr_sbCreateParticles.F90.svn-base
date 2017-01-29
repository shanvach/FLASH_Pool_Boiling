!!****if* source/Grid/GridStructures/withTriangles/pointTopoint/bitmap/gr_sbCreateParticles
!!
!! NAME
!!  gr_sbCreateParticles
!!
!! SYNOPSIS
!!
!!  gr_sbCreateParticles()
!!  
!! DESCRIPTION 
!!  
!!  This routine is called from Grid_initDomain and Driver_evolveFlash for
!!  moving body. The master processor creates 
!!  particles for each triangle in solid body. Also gets the position coordinates in the
!!  particle structure. 
!!
!!
!! ARGUMENTS 
!!
!!***

#include "constants.h"
#include "Flash.h"

Subroutine gr_sbCreateParticles()
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, aelem, totalPart, sumPart
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkBoundBox
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_data, ONLY : gr_meshMe, gr_meshComm

  implicit none
  include "Flash_mpi.h"
  real :: xpos, ypos, zpos
  integer, save :: numParticles
  integer :: b, i, j, p
  logical, parameter :: localLogFile = .true.

  integer :: tmp, numTriangles, ie 
  real :: ei, ej, N1, N2, N3
  integer, parameter :: nSideParticles = 4

  integer :: sumPart_pr,ierr

  call Timers_start("body_create_particles")

  p = 0
  numParticles = (nSideParticles**2 + nSideParticles) / 2
  if (numParticles <= 0) then
     call Driver_abortFlash("Invalid solid body particle distribution")
  end if

  do b = 1, gr_sbNumBodies
     p = 0

     !The solid body master processor creates the particles
     if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then


!        numTriangles = size(gr_sbBodyInfo(b) % triangleCentroids,2)
        numTriangles = gr_sbBodyInfo(b)%NumAelem
        totalPart = numParticles*numTriangles
        gr_sbBodyInfo(b) % totalPart = totalPart

        allocate(gr_sbBodyInfo(b) % particles(1:NPART_PROPS,numParticles*numTriangles))
        gr_sbBodyInfo(b) % particles(:,:) = NONEXISTENT

        !Lattice distribution - copied from pt_initPositionsLattice

        do ie = 1, numTriangles
           do j = 1, nSideParticles
              tmp = 1 - j + nSideParticles
              do i = 1, tmp
                 ei  = real(i-1)/real(nSideparticles-1)
                 ej  = real(j-1)/real(nSideParticles-1)

                    ! Shape functions:
                 N1 = 1 - ei - ej
                 N2 = ei
                 N3 = ej
  
                    ! x,y,z positions of internal particle:
                 xpos = N1*gr_sbBodyInfo(b) % xb(aelem(1,ie)) + N2*gr_sbBodyInfo(b) &
                      % xb(aelem(2,ie)) + N3*gr_sbBodyInfo(b) % xb(aelem(3,ie))
                 ypos = N1*gr_sbBodyInfo(b) % yb(aelem(1,ie)) + N2*gr_sbBodyInfo(b) &
                      % yb(aelem(2,ie)) + N3*gr_sbBodyInfo(b) % yb(aelem(3,ie))
                 zpos = N1*gr_sbBodyInfo(b) % zb(aelem(1,ie)) + N2*gr_sbBodyInfo(b) &
                      % zb(aelem(2,ie)) + N3*gr_sbBodyInfo(b) % zb(aelem(3,ie))

                 p = p + 1
                 
                 gr_sbBodyInfo(b) % particles(BLK_PART_PROP,p) = UNKNOWN
                 gr_sbBodyInfo(b) % particles(POSX_PART_PROP,p) = xpos
                 gr_sbBodyInfo(b) % particles(POSY_PART_PROP,p) = ypos
                 gr_sbBodyInfo(b) % particles(POSZ_PART_PROP,p) = zpos
                 gr_sbBodyInfo(b) % particles(GLOB_PART_PROP,p) = p !local particle counter in solid body.
              end do
           end do
        end do

     else

       gr_sbBodyInfo(b) % totalPart = 0

     end if
  end do

 ! Get Sum of Particles for all bodies on Processor:
  sumPart_pr = 0
  do b = 1, gr_sbNumBodies
     sumPart_pr = sumPart_pr + gr_sbBodyInfo(b) % totalPart
  enddo
  call mpi_allreduce ( sumPart_pr, sumPart, 1, FLASH_INTEGER, &
                       MPI_MAX, gr_meshComm, ierr )


  call Timers_stop("body_create_particles")

End Subroutine gr_sbCreateParticles
