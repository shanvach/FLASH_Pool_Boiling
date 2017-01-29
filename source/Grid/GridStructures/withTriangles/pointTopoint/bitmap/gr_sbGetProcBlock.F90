!!****if* source/Grid/GridStructures/withTriangles/pointTopoint/bitmap/gr_sbGetProcBlock
!!
!! NAME
!!  gr_sbGetProcBlock
!!
!! SYNOPSIS
!!
!!  gr_sbGetProcBlock()
!!  
!! DESCRIPTION 
!!  
!!  This routine is called from Driver_evolveFlash. The master processor gets the
!!  proc and blkID of each particle given its position.
!!
!!
!! ARGUMENTS 
!!
!!***

#include "constants.h"
#include "Flash.h"

Subroutine gr_sbGetProcBlock()
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, totalPart, gr_sbFirstCall
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_interface, ONLY : gr_xyzToBlockLevel
  use Timers_interface, ONLY : Timers_start, Timers_stop
  
#ifdef FLASH_GRID_PARAMESH
  use bittree, only : amr_identify_block
  use Grid_data, ONLY : gr_meshNumProcs,gr_nblockX, gr_nblockY, gr_nblockZ
#else
  use Grid_data, ONLY : gr_axisNumProcs,gr_globalDomain
  use gr_sbData, ONLY : gr_sbIJKtoProc  
#endif
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none
  include "Flash_mpi.h"

  real, dimension(MDIM) :: particlePosn
  integer, dimension(MDIM) :: ijk
  integer :: b, i, proc, blk, lev, lev_save,blkList(MAXBLOCKS),blkCount
  real :: pd_xmin, pd_xmax, pd_ymin, pd_ymax, pd_zmin, pd_zmax, deltax, deltay, deltaz

  real,    dimension(MDIM) :: L_dir
  integer, dimension(MDIM) :: ijkproc
  logical, dimension(MDIM) :: out_domain
  integer :: idir

#ifndef FLASH_GRID_PARAMESH
#ifdef NONFIXEDBLOCKSIZE
  call Driver_abortFlash("gr_sbGetProcBlock : This routine only works in FIXEDBLOCKSIZE mode for UG.")
#else
  ! Compute Domain length in processor number scale:
  L_dir(:) = 1.
  do idir =1,NDIM
     L_dir(idir)=(gr_globalDomain(HIGH,idir)-gr_globalDomain(LOW,idir))/real(gr_axisNumProcs(idir))
  enddo

#endif
#endif

  call Timers_start("body_getProcBlk")
  
  do b = 1, gr_sbNumBodies

     ! IF fixed body and not the first call:
     if ((gr_sbBodyInfo(b)%sbIsFixed .eq. CONSTANT_ONE) .and. (gr_sbFirstCall .eq. CONSTANT_ZERO)) cycle

     if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then
#ifdef FLASH_GRID_PARAMESH
        call RuntimeParameters_get("lrefine_max", lev_save)
        lev = lev_save
        call RuntimeParameters_get("xmin", pd_xmin)
        call RuntimeParameters_get("xmax", pd_xmax)
        !deltax = abs(pd_xmin - pd_xmax)/gr_sbNumBodies
        deltax = 1.5*(pd_xmax-pd_xmin)/real(gr_nblockX*NXB*2**(lev-1))
        if (NDIM >= 2) then
           call RuntimeParameters_get("ymin", pd_ymin)
           call RuntimeParameters_get("ymax", pd_ymax)
           !deltay = abs(pd_ymin - pd_ymax)/gr_sbNumBodies
           deltay = 1.5*(pd_ymax-pd_ymin)/real(gr_nblockY*NYB*2**(lev-1))
        end if
        if(NDIM == 3) then
           call RuntimeParameters_get("zmin", pd_zmin)
           call RuntimeParameters_get("zmax", pd_zmax)
           !deltaz = abs(pd_zmin - pd_zmax)/gr_sbNumBodies
           deltaz = 1.5*(pd_zmax-pd_zmin)/real(gr_nblockZ*NZB*2**(lev-1))
        end if
#endif
           !for each particle get the proc and blk ID using bitmap
        
        totalPart = gr_sbBodyInfo(b) % totalPart
        do i = 1, totalPart
           particlePosn(IAXIS) = gr_sbBodyInfo(b) % particles(POSX_PART_PROP,i)
           if(NDIM >1) then
              particlePosn(JAXIS) = gr_sbBodyInfo(b) % particles(POSY_PART_PROP,i)
           endif
           if(NDIM >2) then
              particlePosn(KAXIS) = gr_sbBodyInfo(b) % particles(POSZ_PART_PROP,i)
           endif
#ifdef FLASH_GRID_PARAMESH
           call gr_xyzToBlockLevel(lev, particlePosn(1:NDIM), ijk(1:NDIM)) 
           call amr_identify_block(gr_meshNumProcs, lev, ijk(1:NDIM), proc, blk)
         
#else

           ! Compute i, j, k positions of processor within proc grid:
           ijkproc(:) = 1
           out_domain(:) = .false.
           do idir =1,NDIM
             ijkproc(idir)=int((particlePosn(idir)-gr_globalDomain(LOW,idir))/L_dir(idir))+1  
            
             if ((ijkproc(idir) .lt. 1) .or. (ijkproc(idir) .gt. gr_axisNumProcs(idir))) out_domain(idir) = .true.

           enddo

           ! If particle outside the boundaries defer treatment:           
           if (any(out_domain)) then
             proc = 0             
             blk  = 0
           else ! Find processor it belongs to.
             proc = gr_sbIJKtoProc(ijkproc(IAXIS),ijkproc(JAXIS),ijkproc(KAXIS))
             blk  = 1
           endif   

#endif
           
           gr_sbBodyInfo(b) % particles(PROC_PART_PROP,i) = proc
           gr_sbBodyInfo(b) % particles(BLK_PART_PROP,i) = blk
           !Check if particle is in the border of the physical domain bound
#ifdef FLASH_GRID_PARAMESH
           if (gr_sbBodyInfo(b) % particles(BLK_PART_PROP,i) < 1) then              

              if (particlePosn(IAXIS) .le. pd_xmin) then
                 particlePosn(IAXIS) = particlePosn(IAXIS) + deltax
              elseif (particlePosn(IAXIS) .ge. pd_xmax) then
                 particlePosn(IAXIS) = particlePosn(IAXIS) - deltax
              end if
              if (NDIM > 1) then
                 if (particlePosn(JAXIS) .le. pd_ymin) then
                    particlePosn(JAXIS) = particlePosn(JAXIS) + deltay
                 elseif (particlePosn(JAXIS) .ge. pd_ymax) then
                    particlePosn(JAXIS) = particlePosn(JAXIS) - deltay
                 end if
              end if
              if (NDIM > 2) then
                 if (particlePosn(KAXIS) .le. pd_zmin) then
                    particlePosn(KAXIS) = particlePosn(KAXIS) + deltaz
                 elseif (particlePosn(KAXIS) .ge. pd_zmax) then
                    particlePosn(KAXIS) = particlePosn(KAXIS) - deltaz
                 end if
              end if
              lev = lev_save ! Reset the Level to lrefine_max, as it might have
                             ! been changed on the previous call to
                             ! amr_identify_block 
              call gr_xyzToBlockLevel(lev, particlePosn(1:NDIM), ijk(1:NDIM))
              call amr_identify_block(gr_meshNumProcs, lev, ijk(1:NDIM), proc, blk)
              gr_sbBodyInfo(b) % particles(PROC_PART_PROP,i) = proc
              gr_sbBodyInfo(b) % particles(BLK_PART_PROP,i) = blk

              if (gr_sbBodyInfo(b) % particles(BLK_PART_PROP,i) < 1) then
                 print *, "less than 1", particlePosn(IAXIS), particlePosn(JAXIS), particlePosn(KAXIS)
              end if
           end if
#endif
        enddo
     endif


  enddo

  call Timers_stop("body_getProcBlk")
     
End Subroutine gr_sbGetProcBlock
