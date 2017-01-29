!!****if* source/Grid/GridStructures/withTriangles/gr_sbCreateParticles
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
!!  particle structure. Finds whether the particle belongs to the master processor.
!!
!!
!! ARGUMENTS 
!!
!!***

#include "constants.h"
#include "Flash.h"

Subroutine gr_sbCreateParticles()
  use Logfile_interface, ONLY : Logfile_open, Logfile_close
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, &
       gr_sbPtNumX, gr_sbPtNumY, gr_sbPtNumZ, gr_sbDebug, &
       aelem, totalPart
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkBoundBox
  use Timers_interface, ONLY : Timers_start, Timers_stop
  implicit none
  include "Flash_mpi.h"
  real, dimension(2,MDIM) :: sb
  real :: xpos, ypos, zpos
  integer, save :: numParticles
  integer :: b, i, j, p
  integer :: logUnit
  logical, parameter :: localLogFile = .true.

  real, dimension(MDIM) :: particlePosn
  real, dimension(2,MDIM) :: boundBox
  integer,dimension(MAXBLOCKS) :: listOfBlocks
  integer :: count, blkID, blkCount
  integer :: tmp, numTriangles, ie 
  real :: ei, ej, N1, N2, N3
  integer, parameter :: nSideParticles = 4

  call Timers_start("body_create_particles")

  p = 0
  numParticles = (nSideParticles**2 + nSideParticles) / 2
  if (numParticles <= 0) then
     call Driver_abortFlash("Invalid solid body particle distribution")
  end if

  do b = 1, gr_sbNumBodies
     p = 0
     if (gr_sbBodyInfo(b) % comm /= MPI_COMM_NULL) then
        numTriangles = size(gr_sbBodyInfo(b) % triangleCentroids,2)
        allocate(gr_sbBodyInfo(b) % particles(1:NPART_PROPS,numParticles*numTriangles))
        totalPart = numParticles*numTriangles
        gr_sbBodyInfo(b) % particles(:,:) = NONEXISTENT

        !The solid body master processor creates the particles
        if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then
           !Lattice distribution - copied from pt_initPositionsLattice

           sb = gr_sbBodyInfo(b) % boundBox

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

                    !For matlab vis.
                    !call Logfile_open(logUnit,localLogFile)
                    !write(logUnit,'(3es14.6)') xpos, ypos, zpos
                    !call Logfile_close(localLogFile)

                    if (gr_sbDebug) then
                       call Logfile_open(logUnit,localLogFile)
                       write(logUnit,'(a,i8,a,3es14.6)') &
                            "Create particle ", p, &
                            " at position ", xpos, ypos, zpos
                       call Logfile_close(localLogFile)
                    end if

                    gr_sbBodyInfo(b) % particles(BLK_PART_PROP,p) = UNKNOWN
                    gr_sbBodyInfo(b) % particles(POSX_PART_PROP,p) = xpos
                    gr_sbBodyInfo(b) % particles(POSY_PART_PROP,p) = ypos
                    gr_sbBodyInfo(b) % particles(POSZ_PART_PROP,p) = zpos
                    gr_sbBodyInfo(b) % particles(GLOB_PART_PROP,p) = p !local particle counter in solid body.
                    call Grid_getListOfBlocks(LEAF, listOfBlocks, count)
                    do blkCount = 1, count
                       blkID = listofBlocks(blkCount)
                       call Grid_getBlkBoundBox(blkID, boundBox)
                       particleposn(IAXIS) = xpos
                       particleposn(JAXIS) = ypos
                       particleposn(KAXIS) = zpos
                       if (all(particleposn(1:NDIM) >= boundBox(LOW,1:NDIM) .and. &
                            particleposn(1:NDIM) < boundBox(HIGH,1:NDIM))) then
                          !particle in Master PE
                          gr_sbBodyInfo(b) % particles(PROC_PART_PROP,p) = &
                               real(gr_sbBodyInfo(b) % myPE)
                          gr_sbBodyInfo(b) % particles(BLK_PART_PROP,p) = blkID
                          exit
                       endif
                    enddo
                    if (blkCount > count) then
                       gr_sbBodyInfo(b) % particles(PROC_PART_PROP,p) = NONEXISTENT
                       gr_sbBodyInfo(b) % particles(BLK_PART_PROP,p) = NONEXISTENT
                    endif
                 end do
              end do
           end do
        end if
     end if
  end do

  call Timers_stop("body_create_particles")

End Subroutine gr_sbCreateParticles
