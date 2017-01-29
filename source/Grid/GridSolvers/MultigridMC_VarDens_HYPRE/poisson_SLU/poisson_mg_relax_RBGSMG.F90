!*******************************************************************************

!  Routine:     mg_relax()

!  Description: Perform a number of smoothing iterations for multigrid.
!               This version implements Gauss-Seidel with red-black ordering 
!               using the Poisson equation.

!  Parameters:  level       Level to relax on.
!               irhs        Right-hand side (source) of equation.
!               ilhs        Left-hand side (solution) of equation.  Receives
!                           the smoothed result.
!               nsmooth     Number of iterations to apply.


subroutine poisson_mg_relax_RBGSMG (level, irhs, ilhs, nsmooth, idenvar)

  !===============================================================================


  use mg_common, ONLY:  ili, iui, jli, jui, kli, kui, gr_mgDiffOpDiscretize

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getBlkPtr,       &
                                Grid_releaseBlkPtr,   &
                                Grid_getBlkBoundBox,  &
                                Grid_getBlkIndexLimits,&
                                Grid_getListOfBlocks

  use Grid_data, ONLY : gr_meshMe

  use workspace, ONLY: work
  use tree, only : nodetype,lrefine
  use paramesh_dimensions, only : nguard_work

  implicit none
#include "Multigrid.h"
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer       :: level, irhs, ilhs, nsmooth, idebug, idenvar
  integer       :: firstpass = 1

  integer       :: iter, i, j, k, lb, ierr, Neff, lnblocks
  integer       :: redblackpass, isweep, jsweep, ksweep
  real          :: lerrnorm(2), errnorm(2), nxbinv, nybinv, nzbinv, rTest, rTest1
  real          :: delx, dely, delz, error, prefac, critfac, pi
  real          :: idelx,idely,idelz

  logical       :: done

  character(len=256) :: str_buffer

  integer, parameter :: iterating_to_convergence_limit = 0
  integer, parameter :: MAXDIMS = 3
  real, dimension(MAXDIMS) :: size

  logical, save :: first_call = .true.
  real, save    :: mgrid_smooth_tol
  integer, save :: Nblockx, Nblocky, Nblockz, myPE, masterPE

  real, pointer, dimension(:,:,:,:), save :: unk

  !- kpd -
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData, &
                                       facevarx, facevary, facevarz
  !- kpd - 
  real :: MdensXL, MdensXR, MdensYL, MdensYR, MdensZL, MdensZR


  integer, parameter :: nxb = NXB
  integer, parameter :: nyb = NYB
  integer, parameter :: nzb = NZB
  integer, parameter :: ndim = NDIM


  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox

  integer blockcount,ii,jj,blockID

  real bsize(MDIM),coord(MDIM)

  integer blockList(MAXBLOCKS)

  !               On the first call, check to make sure that the requested
  !               convergence criterion is reasonable.  The smoother relaxes
  !               long wavelengths (relative to the grid) most slowly.  The
  !               termination criterion stops the iteration when the iterate
  !               has stopped changing to within some tolerance.  This
  !               tolerance should not be larger than the convergence rate
  !               for the longest-wavelength mode.  The convergence rate
  !               used here is given by Briggs et al. (_A Multigrid Tutorial_)
  !               for Dirichlet boundary conditions.

  !               The way we do this comparison at the end may not be quite
  !               right... should probably compute the ratio of the new
  !               residual norm to the old residual norm.

!  call timer_start("relax")

  error    = 1.

  if (first_call) then
     call RuntimeParameters_get('mgrid_smooth_tol',    mgrid_smooth_tol)
     call RuntimeParameters_get('Nblockx',  Nblockx)
     call RuntimeParameters_get('Nblocky',  Nblocky)
     call RuntimeParameters_get('Nblockz',  Nblockz)
     myPE     = gr_meshMe 
     masterPE = MASTER_PE 
  end if

  if (first_call) then
     if (MyPE == MasterPE) then
        pi = PI
        i = 1      ! Will only try to converge via relaxation on level 1
        Neff = max( nxb*Nblockx*2**(i-1), & 
             &                  nyb*Nblocky*2**(i-1), & 
             &                  nzb*Nblockz*2**(i-1) )
        critfac = 2.*0.66667*sin(pi/(2.*Neff))**2
        if (mgrid_smooth_tol >= critfac) then
           write (str_buffer,*) 'mg_relax: termination', & 
                & ' tolerance of ', mgrid_smooth_tol, ' is larger', &
                & ' than the slowest', & 
                & ' convergence rate on level ', i, ' (', critfac, ').', &
                & ' Long-wavelength modes', & 
                & ' may not be able to converge on this level.'
           write(*,*) str_buffer
           !call stamp_logfile(str_buffer, 'warning')
        endif
     endif
     first_call = .false.
  endif

  !               Iteration loop.

  call Grid_getLocalNumBlks(lnblocks)

  iter    = 0
  nxbinv  = 1./nxb**2
  nybinv  = 1./nyb**2
  nzbinv  = 1./nzb**2
  done    = .false.

  
  !               Update boundary zones of the LHS array for this iteration.
  !               Currently the mesh package doesn't supply a general boundary
  !               update routine, so we must copy the LHS into the "work" array,
  !               then update its boundaries.  This is handled by mg_bndry().
  !               This relaxer assumes the LHS variable is defined on all blocks
  !               (not just leaf nodes), so we call mg_bndry with leaf_only == 0.
     
  ! copy unk into work and fill work's guardcells
! ***

  call gr_mgBndry (level, ilhs, 2, 0, MG_COPY_UNK_TO_WORK, MG_BEGIN_SERIES) 

  do while ( .not. done )

     if (iter > iterating_to_convergence_limit) then 
        ! if we're iterating to convergence, 
        ! copy the old iterate for convergence testing.  
        do lb = 1, lnblocks
           if (lrefine(lb) == level) then

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unk,CENTER)

              do k = kli, kui
                 do j = jli, jui
                    do i = ili, iui
                       unk(ilhs,i,j,k) = work(i,j,k,lb,1)
                    enddo
                 enddo
              enddo

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unk,CENTER)


           endif
        enddo
     end if

     !               Perform the Gauss-Seidel iteration step. 
 
!     call timer_start("relax calc")

     if (ndim == 1) then

     elseif (ndim == 2) then

        do lb = 1, lnblocks
           if (lrefine(lb) == level) then
              call Grid_getBlkPhysicalSize(lb,size)

              !- kpd - 
              call Grid_getBlkPtr(lb,facexData,FACEX)
              call Grid_getBlkPtr(lb,faceyData,FACEY)

              !- kpd - These are all 1/dx^2 or 1/dy^2
              delx = size(1)**2 * nxbinv
              dely = size(2)**2 * nybinv
              prefac = 0.5 / (delx+dely)

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unk,CENTER)
              jsweep = 0
              do redblackpass = 1, 2
                 isweep = jsweep
                 do j = jli, jui
                    do i = ili + isweep, iui, 2

                       !-------------------------------------------------------
                       !- kpd - 2nd Order Central Poisson Solve in X- and Y-dir
                       !-------------------------------------------------------

                       !- kpd - Added for variable density. Calculate face densities...
                       MdensXL = facevarx(idenvar,i  ,j  ,1)
                       MdensXR = facevarx(idenvar,i+1,j  ,1)
                       MdensYL = facevary(idenvar,i  ,j  ,1)
                       MdensYR = facevary(idenvar,i  ,j+1,1)

                       if (MdensXL==0. .OR. MdensXR==0. .OR. MdensYL==0. .OR. MdensYR==0. ) then
                          print*,"ERROR poisson_mg_relax_RBGSMG.F90: Density Equals Zero at a FACE",lb,i,j,MdensXR,MdensYR
                       elseif (MdensXL.gt.1000.0 .OR. MdensXR.gt.1000.0 .OR. MdensYL.gt.1000.0 .OR. MdensYR.gt.1000.0 ) then
                          print*,"ERROR poisson_mg_relax_RBGSMG.F90: Inverse Density Greater than 1.0",lb,i,j,MdensXR,MdensYR
                       end if


                       !- kpd - Added for variable density...
                       !        div (1/rho) * (grad(P)) = dx^2 * (div(uStar + sig*K))
                       !        Note: unk(irhs,...) is the residual (r=RHS-A*P)
                       !        Note: work is the error correction (Ae=r, Pnew = Pold + e)
                       !        Note: delx is dx^2
                       work(i,j,kli:kui,lb,1) = &
                          ( MdensXR/delx * work(i+1,j,kli:kui,lb,1) +   &
                            MdensXL/delx * work(i-1,j,kli:kui,lb,1) +   &
                            MdensYR/dely * work(i,j+1,kli:kui,lb,1) +   &
                            MdensYL/dely * work(i,j-1,kli:kui,lb,1) -   &
                            unk(irhs,i,j,kli:kui)             ) / &
                          ( MdensXR/delx + MdensXL/delx       +   &
                            MdensYR/dely + MdensYL/dely       )

                    enddo
                    isweep = 1 - isweep
                 enddo
                 jsweep = 1 - jsweep
              enddo

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unk,CENTER)

              call Grid_releaseBlkPtr(lb,facexData,FACEX)
              call Grid_releaseBlkPtr(lb,faceyData,FACEY)

           endif
        enddo

     else

     endif

     !               Now compute the error norm (change in LHS from one iteration
     !               to the next).  This will be used to limit the number of
     !               iterations (see above).  Here we use the maximum norm rather
     !               than the L2 norm computed by mg_norm.

     !               However, don't start checking the convergence criterion
     !               until we've done ten iterations.  This saves work when doing
     !               the regular smoothing on fine levels, which isn't expected
     !               to produce convergence by itself.


     iter = iter + 1

     if (iter > iterating_to_convergence_limit) then

        lerrnorm(:) = 0.
      
        do lb = 1, lnblocks
           if (lrefine(lb) == level) then

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unk,CENTER)

              do k = kli, kui
                 do j = jli, jui
                    do i = ili, iui
                       lerrnorm(1) = & 
                            max( lerrnorm(1), & 
                            abs(unk(ilhs,i,j,k)-work(i,j,k,lb,1)) )
                       lerrnorm(2) = & 
                            max( lerrnorm(2), abs(unk(ilhs,i,j,k)) )
                    enddo
                 enddo
              enddo

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unk,CENTER)

           endif
        enddo

        call mpi_allreduce (lerrnorm, errnorm, 2, MPI_DOUBLE_PRECISION, & 
             MPI_MAX, MPI_COMM_WORLD, ierr)

        error = mgrid_smooth_tol*errnorm(2) - errnorm(1)

     endif

     !               End of loop.  Iterate until nsmooth iterations have been
     !               performed, or until the convergence criterion has been met.

     done = (iter == nsmooth) .or. ((iter > iterating_to_convergence_limit) .and. (error > 0.))
!     call timer_stop("relax calc")


!     if (done .and. (gr_meshMe .eq. 0)) write(*,*) 'Relax Iterations=',iter 

     ! fill work's guardcells
     if (.not. done) &
     call gr_mgBndry (level, ilhs, nguard_work, 0, MG_EXCHANGE_WORK, MG_CONTINUE_SERIES) 

  enddo


! ***
! Why is this needed ?  It is called above !!!
!!!  call mg_bndry(level, ilhs, nguard_work, 0, EXCHANGE_WORK, CONTINUE_SERIES)
  call gr_mgBndry(level, ilhs, nguard_work, 0, MG_EXCHANGE_WORK, MG_END_SERIES)

  do lb = 1, lnblocks
     if (lrefine(lb) == level) then

        ! Point to blocks center vars:
        call Grid_getBlkPtr(lb,unk,CENTER)

        !print*,"KPD LB",lb

        do k = kli, kui
           do j = jli, jui
              do i = ili, iui
                 unk(ilhs,i,j,k) = work(i,j,k,lb,1)
              enddo
           enddo
        enddo

        ! Point to blocks center vars:
        call Grid_releaseBlkPtr(lb,unk,CENTER)


     endif
  enddo

  !===============================================================================

!  call timer_stop("relax")
  return
end subroutine poisson_mg_relax_RBGSMG
