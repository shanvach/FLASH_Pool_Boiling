!!****if* source/Simulation/SimulationMain/INavierStokes/2D/IB_Cyl_parallelIBVP/old/gr_sbCreateParticles
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
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, aelem, totalPart,sumPart
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkBoundBox
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_data, ONLY : gr_meshMe, gr_meshComm

  use ib_interface, ONLY : ib_countParticles, ib_mapParticles

  implicit none
  include "Flash_mpi.h"
  real :: xpos, ypos, zpos
  integer, save :: numParticles
  integer :: b, i, j, p
  logical, parameter :: localLogFile = .true.

  integer :: sumPart_pr,ierr

  integer :: tmp, numTriangles, ie 
  real :: ei, ej, N1, N2, N3
  integer, parameter :: nSideParticles = 2

  real :: xvel,yvel,zvel,area,areai
  integer :: ael_1,ael_2

  
  call Timers_start("body_create_particles")

!!$  p = 0
!!$#if NDIM == 2
!!$  numParticles = nSideParticles 
!!$#elif NDIM == 3
!!$  numParticles = (nSideParticles**2 + nSideParticles) / 2
!!$#endif
!!$  if (numParticles <= 0) then
!!$     call Driver_abortFlash("Invalid solid body particle distribution")
!!$  end if

  do b = 1, gr_sbNumBodies

!!$     p = 0
!!$     numTriangles = gr_sbBodyInfo(b)%NumAelem  

!     size(gr_sbBodyInfo(b) % triangleCentroids,2)
!     allocate(gr_sbBodyInfo(b) % particles(1:NPART_PROPS,numParticles*numTriangles))
!     totalPart = numParticles*numTriangles!3810
!     gr_sbBodyInfo(b) % totalPart = totalPart 
!     gr_sbBodyInfo(b) % particles(:,:) = NONEXISTENT

     !The solid body master processor creates the particles
     if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then


        ! Count the number of particles necessary for this body:
        call ib_countParticles(b,totalPart)
        gr_sbBodyInfo(b) % totalPart = totalPart

        ! Allocate the particles array:
        allocate(gr_sbBodyInfo(b) % particles(1:NPART_PROPS,gr_sbBodyInfo(b)%totalPart))
        gr_sbBodyInfo(b) % particles(:,:) = NONEXISTENT


        ! Now do the Particle Mapping:
        call ib_mapParticles(b)
 
!!$        do ie = 1, numTriangles
!!$#if NDIM == 2
!!$
!!$           ael_1 = gr_sbBodyInfo(b) % AELEM(2,ie) !aelem(1,ie)
!!$           ael_2 = gr_sbBodyInfo(b) % AELEM(3,ie) !aelem(2,ie)
!!$
!!$           area = sqrt( (gr_sbBodyInfo(b) % xb(ael_2) - gr_sbBodyInfo(b) % xb(ael_1))**2. + &
!!$                        (gr_sbBodyInfo(b) % yb(ael_2) - gr_sbBodyInfo(b) % yb(ael_1))**2.)
!!$
!!$           do i = 1, nSideParticles
!!$
!!$              ei  = real(i-1)/real(nSideparticles-1)
!!$
!!$              ! Shape functions:
!!$              N1 = 1 - ei
!!$              N2 = ei
!!$                 
!!$              ! x,y,z positions of internal particle:
!!$              xpos = N1*gr_sbBodyInfo(b) % xb(ael_1) + N2*gr_sbBodyInfo(b) &
!!$                   % xb(ael_2) 
!!$              ypos = N1*gr_sbBodyInfo(b) % yb(ael_1) + N2*gr_sbBodyInfo(b) &
!!$                      % yb(ael_2)
!!$              zpos = 0.
!!$
!!$              !write(*,*) xpos,ypos,gr_sbBodyInfo(b)%yb(ael_1),N2*gr_sbBodyInfo(b)%yb(ael_2),ael_1,ael_2
!!$
!!$              ! x,y,z velocities of internal particles
!!$              xvel =  N1*gr_sbBodyInfo(b) % ubd(ael_1) + N2*gr_sbBodyInfo(b) &
!!$                      % ubd(ael_2) 
!!$              yvel =  N1*gr_sbBodyInfo(b) % vbd(ael_1) + N2*gr_sbBodyInfo(b) &
!!$                      % vbd(ael_2)
!!$              zvel = 0.
!!$
!!$
!!$              ! Area associated with particle:
!!$              areai = area/real(nSideParticles-1)
!!$              if ((i .eq. 1).or.(i .eq. nSideParticles)) areai = areai/2.
!!$              
!!$              ! Dump into Particles:
!!$              p = p + 1
!!$                 
!!$              gr_sbBodyInfo(b) % particles(BLK_PART_PROP,p) = UNKNOWN
!!$              gr_sbBodyInfo(b) % particles(POSX_PART_PROP,p) = xpos
!!$              gr_sbBodyInfo(b) % particles(POSY_PART_PROP,p) = ypos
!!$              gr_sbBodyInfo(b) % particles(POSZ_PART_PROP,p) = zpos
!!$              gr_sbBodyInfo(b) % particles(VELX_PART_PROP,p) = xvel
!!$              gr_sbBodyInfo(b) % particles(VELY_PART_PROP,p) = yvel
!!$              gr_sbBodyInfo(b) % particles(VELZ_PART_PROP,p) = zvel
!!$              gr_sbBodyInfo(b) % particles(FUL_PART_PROP,p) = 0.
!!$              gr_sbBodyInfo(b) % particles(FVL_PART_PROP,p) = 0. 
!!$              gr_sbBodyInfo(b) % particles(FWL_PART_PROP,p) = 0.
!!$              gr_sbBodyInfo(b) % particles(TAG_PART_PROP,p) = 1.
!!$              gr_sbBodyInfo(b) % particles(AREA_PART_PROP,p) = areai
!!$              gr_sbBodyInfo(b) % particles(GLOB_PART_PROP,p) = p !local particle counter in solid body.
!!$           end do
!!$#elif NDIM == 3
!!$           do j = 1, nSideParticles
!!$              tmp = 1 - j + nSideParticles
!!$              do i = 1, tmp
!!$                 ei  = real(i-1)/real(nSideparticles-1)
!!$                 ej  = real(j-1)/real(nSideParticles-1)
!!$
!!$                    ! Shape functions:
!!$                 N1 = 1 - ei - ej
!!$                 N2 = ei
!!$                 N3 = ej
!!$  
!!$                    ! x,y,z positions of internal particle:
!!$                 xpos = N1*gr_sbBodyInfo(b) % xb(aelem(1,ie)) + N2*gr_sbBodyInfo(b) &
!!$                      % xb(aelem(2,ie)) + N3*gr_sbBodyInfo(b) % xb(aelem(3,ie))
!!$                 ypos = N1*gr_sbBodyInfo(b) % yb(aelem(1,ie)) + N2*gr_sbBodyInfo(b) &
!!$                      % yb(aelem(2,ie)) + N3*gr_sbBodyInfo(b) % yb(aelem(3,ie))
!!$                 zpos = N1*gr_sbBodyInfo(b) % zb(aelem(1,ie)) + N2*gr_sbBodyInfo(b) &
!!$                      % zb(aelem(2,ie)) + N3*gr_sbBodyInfo(b) % zb(aelem(3,ie))
!!$
!!$                 p = p + 1
!!$                 
!!$                 gr_sbBodyInfo(b) % particles(BLK_PART_PROP,p) = UNKNOWN
!!$                 gr_sbBodyInfo(b) % particles(POSX_PART_PROP,p) = xpos
!!$                 gr_sbBodyInfo(b) % particles(POSY_PART_PROP,p) = ypos
!!$                 gr_sbBodyInfo(b) % particles(POSZ_PART_PROP,p) = zpos
!!$                 gr_sbBodyInfo(b) % particles(GLOB_PART_PROP,p) = p !local particle counter in solid body.
!!$              end do
!!$           end do
!!$#endif
!!$        end do

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
