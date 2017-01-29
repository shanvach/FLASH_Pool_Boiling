!!****if* source/Grid/GridStructures/gr_sbCreateParticles
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
!!  This routine is called from Grid_initDomain. The master processor creates 
!!  particles for each solid body. Also gets the position coordinates in the
!!  particle structure.
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
       gr_sbPtNumX, gr_sbPtNumY, gr_sbPtNumZ, gr_sbDebug
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkBoundBox
  use Timers_interface, ONLY : Timers_start, Timers_stop
  implicit none
  include "Flash_mpi.h"
  real, dimension(2,MDIM) :: sb
  real, dimension(MDIM) :: dxParticle
  real :: xpos, ypos, zpos
  integer, save :: numParticles
  integer :: b, i, j, k, p
  integer :: logUnit
  logical, parameter :: localLogFile = .true.

  real, dimension(MDIM) :: particlePosn
  real, dimension(2,MDIM) :: boundBox
  integer,dimension(MAXBLOCKS) :: listOfBlocks
  integer :: count, blkID, blkCount

  call Timers_start("body_create_particles")

  p = 0
  numParticles = gr_sbPtNumX
  if (NDIM >= 2) numParticles = numParticles * gr_sbPtNumY
  if (NDIM == 3) numParticles = numParticles * gr_sbPtNumZ
  if (numParticles <= 0) then
     call Driver_abortFlash("Invalid solid body particle distribution")
  end if

  do b = 1, gr_sbNumBodies
     p = 0
     if (gr_sbBodyInfo(b) % comm /= MPI_COMM_NULL) then
        allocate(gr_sbBodyInfo(b) % particles(1:NPART_PROPS,numParticles))
        gr_sbBodyInfo(b) % particles(:,:) = NONEXISTENT

        !The solid body master processor creates the particles
        if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then
           !Lattice distribution - copied from pt_initPositionsLattice

           sb = gr_sbBodyInfo(b) % boundBox

           ! determine particle spacing as dxParticle
           dxParticle = 0.0   ! initialize to zero for all dimensions
           dxParticle(IAXIS) = (sb(HIGH,IAXIS) - sb(LOW,IAXIS)) / gr_sbPtNumX
           if (NDIM >= 2) then
              dxParticle(JAXIS) = (sb(HIGH,JAXIS) - sb(LOW,JAXIS)) / gr_sbPtNumY
           end if
           if (NDIM == 3) then
              dxParticle(KAXIS) = (sb(HIGH,KAXIS) - sb(LOW,KAXIS)) / gr_sbPtNumZ
           end if

           !! initialization in case of lower dimensionality
           zpos = 0.0
           ypos = 0.0
           xpos = 0.0

           loop_x:  do i = 1, gr_sbPtNumX
              xpos = (i-0.5)*dxParticle(IAXIS) + sb(LOW,IAXIS)

              loop_y:  do j = 1, gr_sbPtNumY
                 if (NDIM >= 2) then
                    ypos = (j-0.5)*dxParticle(JAXIS) + sb(LOW,JAXIS)
                 end if

                 loop_z:  do k = 1, gr_sbPtNumZ
                    if (NDIM == 3) then
                       zpos = (k-0.5)*dxParticle(KAXIS) + sb(LOW,KAXIS)
                    end if

                    p = p + 1

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
                 end do loop_z
              end do loop_y
           end do loop_x
        end if
     end if
  end do

  call Timers_stop("body_create_particles")

End Subroutine gr_sbCreateParticles
