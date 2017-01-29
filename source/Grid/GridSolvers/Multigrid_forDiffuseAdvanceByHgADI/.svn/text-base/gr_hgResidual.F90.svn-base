!!****if* source/Grid/GridSolvers/Multigrid_forDiffuseAdvanceByHgADI/gr_hgResidual
!!
!! NAME
!!  gr_hgResidual
!!
!! SYNOPSIS
!!  
!!  gr_hgResidual(integer, intent(in) :: level
!!                integer, intent(in) :: gr_iSource
!!                integer, intent(in) :: gr_iSoln
!!                integer, intent(in) :: ires)
!!
!! DESCRIPTION
!!  
!!  Compute the current residual for the solution to the Poisson equation.
!!  This routine uses the work array and boundary filling on the work array
!!  as the boundary values are used in the computation of the residual one
!!  block in.  Only the leaf block residuals are calculated.
!!
!!  This is done by applying the first-order laplacian stencil as defined 
!!  by hg_cx, hg_cy, and hg_cz in gr_hgData and gr_hgInit.
!!
!! ARGUMENTS
!!
!!  level        - the level (of leaf blocks) to take the residual for
!!  gr_iSource - the density variable
!!  gr_iSoln   - the potential variable
!!  ires         - the residual variable
!!
!! RESULT
!!
!!  The residual between the source and solutions at the blocks at level 
!!  is placed the variable ires.
!!
!! NOTES
!!
!!  The original description was somewhat dishonest in that extrapolations
!!  are done for the exterior, however the work boundary fill is quite real.
!!  A potential optimization would be to do this ONCE on all leaf blocks,
!!  saving us a number of boundary fills.
!!
!!  Note that gr_iSource is copied into work on the leaf blocks here.
!!
!!***

!!REORDER(5): unk

subroutine gr_hgResidual(level, gr_iSource, gr_iSoln, ires,dt,chi,theta)

!==============================================================================
#include "Flash.h"
#include "Multigrid.h"
#include "constants.h"

  use tree, ONLY : lnblocks,lrefine_min,lrefine,bsize
  use Grid_interface, ONLY : Grid_getDeltas
  use physicaldata, ONLY : unk
  use workspace, ONLY : work
  use gr_hgData, ONLY : hg_ili, hg_iui, hg_jli, hg_jui, hg_kli, hg_kui, &
                        gr_hgBndTypes, gr_hgSaveNodetype

  use Timers_interface, ONLY : Timers_start, Timers_stop
  use gr_hgInterface, ONLY: gr_hgBndry

  implicit none

  include 'Flash_mpi.h'

  integer, intent(in)          :: level, gr_iSource, gr_iSoln, ires
  integer                      :: b, i, j, k, ii, jj, kk, ierr

  real, dimension(MDIM)        :: deltas
  real                         :: K_D
  real                         :: dxinv2, dyinv2, dzinv2
  real,intent(IN),OPTIONAL     :: dt, chi, theta
  real                         :: cond_L,cond_R, LD, MD, UD
  
  
  !=======================================================================

  call Timers_start("gr_hgResidual")
  

  
  !call mg_bndry(level, gr_iSoln, 1, 0, COPY_UNK_TO_WORK, STANDALONE)
  ! Use (EXCHANGE_WORK, CONTINUE_SERIES) under the assumption that we're
  ! calling gr_hgResidual immediately after a call to hg_solveLevel --
  ! performance savings by not filling work() again
  call gr_hgBndry(level, gr_iSoln, NGUARD, 0, &
                MG_COPY_UNK_TO_WORK, MG_STANDALONE, .false.)  !haven't checked this
  
  do b = 1, lnblocks
     if (lrefine(b) == level) then

        call Grid_getDeltas(b,deltas)
                
        do k = NGUARD*K3D+1, NGUARD*K3D+NZB           ! working on interior only
           do j = NGUARD*K2D+1, NGUARD*K2D+NYB
              do i = NGUARD+1, NGUARD+NXB

                 unk(ires,i,j,k,b) = unk(gr_iSource, i,j,k,b)

                 Cond_R = chi ! TO DO : Use COND_VAR
                 Cond_L = chi ! TO DO : Use COND_VAR

                 LD = -theta*Cond_L*(dt/(deltas(1)**2))
                 UD = -theta*Cond_R*(dt/(deltas(1)**2))
                 MD =  1.0 + theta*(Cond_R + Cond_L)*(dt/(deltas(1)**2))

                 unk(ires,i,j,k,b) = unk(ires,i,j,k,b) - (LD*work(i-1,j,k,b,1) + UD*work(i+1,j,k,b,1))

                 if (NDIM .ge. 2) then

                     Cond_R = chi ! TO DO : Use COND_VAR
                     Cond_L = chi ! TO DO : Use COND_VAR

                     MD = MD + theta*(Cond_R + Cond_L)*(dt/(deltas(2)**2))

                     LD = -theta*Cond_L*(dt/(deltas(2)**2))
                     UD = -theta*Cond_R*(dt/(deltas(2)**2))

                     unk(ires,i,j,k,b) = unk(ires,i,j,k,b) - (LD*work(i,j-1,k,b,1) + UD*work(i,j+1,k,b,1))

                 end if

                 if (NDIM .ge. 3) then
 
                     Cond_R = chi ! TO DO : Use COND_VAR
                     Cond_L = chi ! TO DO : Use COND_VAR

                     MD = MD + theta*(Cond_R + Cond_L)*(dt/(deltas(3)**2))

                     LD = -theta*Cond_L*(dt/(deltas(3)**2))
                     UD = -theta*Cond_R*(dt/(deltas(3)**2))

                     unk(ires,i,j,k,b) = unk(ires,i,j,k,b) - (LD*work(i,j,k-1,b,1) + UD*work(i,j,k+1,b,1))
                endif

                unk(ires,i,j,k,b) = unk(ires,i,j,k,b) - (MD*work(i,j,k,b,1))




              enddo
           enddo
        enddo
     end if 
  enddo


  call Timers_stop("gr_hgResidual")

  !=======================================================================
  
  return
end subroutine gr_hgResidual
