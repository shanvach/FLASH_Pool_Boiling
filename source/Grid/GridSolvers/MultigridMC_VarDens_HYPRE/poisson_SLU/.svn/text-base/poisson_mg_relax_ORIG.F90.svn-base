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


subroutine poisson_mg_relax (level, irhs, ilhs, nsmooth)

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

  integer       :: level, irhs, ilhs, nsmooth, idebug
  integer       :: firstpass = 1

  integer       :: iter, i, j, k, lb, ierr, Neff, lnblocks
  integer       :: redblackpass, isweep, jsweep, ksweep
  real          :: lerrnorm(2), errnorm(2), nxbinv, nybinv, nzbinv
  real          :: delx, dely, delz, error, prefac, critfac, pi
  real          :: idelx,idely,idelz
  real          :: MdensXL, MdensXR, MdensYL, MdensYR, MdensZL, MdensZR

  logical       :: done

  character(len=256) :: str_buffer

  integer, parameter :: iterating_to_convergence_limit = 0
  integer, parameter :: MAXDIMS = 3
  real, dimension(MAXDIMS) :: size

  logical, save :: first_call = .true.
  real, save    :: mgrid_smooth_tol
  integer, save :: Nblockx, Nblocky, Nblockz, myPE, masterPE

  real, pointer, dimension(:,:,:,:), save :: unk,facexData,faceyData,facezData

  integer, parameter :: nxb = NXB
  integer, parameter :: nyb = NYB
  integer, parameter :: nzb = NZB
  integer, parameter :: ndim = NDIM


  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox

  integer blockcount,ii,jj

  real bsize(MDIM),coord(MDIM)


  !- kpd - New Local Data
  integer pIndex,iNonZero
  real :: ccAmat_2D(GRID_ILO_GC:GRID_IHI_GC*GRID_JHI_GC,GRID_JLO_GC:GRID_IHI_GC*GRID_JHI_GC)
  real :: ccSol_2D (GRID_ILO_GC:GRID_IHI_GC*GRID_JHI_GC,GRID_JLO_GC:GRID_IHI_GC*GRID_JHI_GC)
  real :: ccRHS_2D (GRID_ILO_GC:GRID_IHI_GC*GRID_JHI_GC)
 

  !==========================================================================

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

!================================================================================================
!================================================================================================
!             - kpd - Iteration loop for the pressure poisson at a given level
!================================================================================================
!================================================================================================

  call Grid_getLocalNumBlks(lnblocks)

  iter    = 0
  nxbinv  = 1./nxb**2
  nybinv  = 1./nyb**2
  nzbinv  = 1./nzb**2
  done    = .false.
  
  ! Update boundary zones of the LHS array for this iteration.
  ! Currently the mesh package doesn't supply a general boundary
  ! update routine, so we must copy the LHS into the "work" array,
  ! then update its boundaries.  This is handled by mg_bndry().
  ! This relaxer assumes the LHS variable is defined on all blocks
  ! (not just leaf nodes), so we call mg_bndry with leaf_only == 0.
  !- kpd - Copy unk into work and fill work's guardcells
  call gr_mgBndry (level, ilhs, 2, 0, MG_COPY_UNK_TO_WORK, MG_BEGIN_SERIES) 

 !=======================
  do while ( .not. done )
 !=======================

     !- kpd - The previous iteration solution (work) is stored 
     !           in ilhs. This is done for convergence testing
     if (iter > iterating_to_convergence_limit) then 
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

     !----------------------------------------- 
     ! Perform the Gauss-Seidel iteration step. 
     !----------------------------------------- 

     ! call timer_start("relax calc")

     !- kpd - 1-D Gauss Seidel
     if (ndim == 1) then
        do lb = 1, lnblocks
           if (lrefine(lb) == level) then
              call Grid_getBlkPhysicalSize(lb,size)
              delx = size(1)**2 * nxbinv

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unk,CENTER)

              isweep = 0
              do redblackpass = 1, 2
                 do i = ili+isweep, iui, 2  
          
                    !- kpd - 2nd Order Central Poisson Solve in X-dir
                    !- kpd - P(i) = 1/2 ( P(i-1) + P(i+1) - dx^2*RHS )
                    work(i,jli:jui,kli:kui,lb,1) = & 
                         &  0.5 * ( work(i-1,jli:jui,kli:kui,lb,1) + & 
                         &          work(i+1,jli:jui,kli:kui,lb,1) - & 
                         &          delx*unk(irhs,i,jli:jui,kli:kui) )
                 end do
                 isweep = 1 - isweep 
              end do

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unk,CENTER)

           endif
        enddo

     !************************************************************************
     !************************************************************************
     !************************************************************************

     !----------------------------------
     !- kpd - 2-D UMF-Pack Direct Solver
     !----------------------------------
     elseif (ndim == 2) then

        do lb = 1, lnblocks

           if (lrefine(lb) == level) then

              call Grid_getBlkPhysicalSize(lb,size)

              !- kpd - For UG case, these are all 1/dx^2
              delx = size(1)**2 * nxbinv
              dely = size(2)**2 * nybinv
              prefac = 0.5 / (delx+dely)

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unk,CENTER)
              !- kpd- Added for variable density formulation
              call Grid_getBlkPtr(lb,facexData,FACEX)
              call Grid_getBlkPtr(lb,faceyData,FACEY)

              !- kpd - Initialize (zero) coeff matrix and RHS before assignment
              pIndex   = 0
              iNonZero = 0
              do j = GRID_JLO_GC,GRID_JHI_GC
                 do i = GRID_ILO_GC,GRID_IHI_GC
                    pIndex = pIndex + 1
                    ccAmat_2D(pIndex,GRID_ILO_GC:GRID_IHI_GC*GRID_JHI_GC) = 0.d0
                    ccAmat_2D(pIndex,pIndex) = 1.d0
                    ccRHS_2D (pIndex) = 0.d0
                    ccSol_2D (i,j) = 0.d0
                 end do
              end do

              !- kpd - Assign coefficient matrix and RHS for directo solver
              do j = jli, jui
                 do i = ili, iui

                    pIndex = i+(j-1)*(GRID_IHI_GC)

                    !- kpd - Added for variable density. Calculate face densities...
                    MdensXL = ( facexData(RH1F_FACE_VAR,i  ,j  ,1) + &
                                facexData(RH2F_FACE_VAR,i  ,j  ,1) )   ! Inverse density on left face.
                    MdensXR = ( facexData(RH1F_FACE_VAR,i+1,j  ,1) + &
                                facexData(RH2F_FACE_VAR,i+1,j  ,1) )   ! Inverse density on right face.
                    MdensYL = ( faceyData(RH1F_FACE_VAR,i  ,j  ,1) + &
                                faceyData(RH2F_FACE_VAR,i  ,j  ,1) )   ! Inverse density on bottom face.
                    MdensYR = ( faceyData(RH1F_FACE_VAR,i  ,j+1,1) + &
                                faceyData(RH2F_FACE_VAR,i  ,j+1,1) )   ! Inverse density on top face.

                    ccAmat_2D(pIndex,pIndex-GRID_IHI_GC) =   MdensYL/dely 
                    ccAmat_2D(pIndex,pIndex-1          ) =   MdensXL/delx 
                    ccAmat_2D(pIndex,pIndex            ) = - MdensXL/delx - MdensXR/delx &
                                                           - MdensYL/dely - MdensYR/dely 
                    ccAmat_2D(pIndex,pIndex+1          ) =   MdensXR/delx
                    ccAmat_2D(pIndex,pIndex+GRID_IHI_GC) =   MdensYR/dely

                    iNonZero = iNonZero + 5

                    ccRHS_2D (pIndex                   ) =   unk(irhs,i,j,1)

                 enddo
              enddo


              !- kpd - Set Boundary Conditions for coeff matrix and RHS 
              do i = ili, iui
                 !- kpd - Bottom cells left to right i
                 ccAmat_2D( (jli-2)*GRID_IHI_GC + i , &
                            (jli-2)*GRID_IHI_GC + i ) = -1.d0

                 !- kpd - Periodic BC
                 !!ccAmat_2D( (jli-2)*GRID_IHI_GC + i + GRID_IHI_GC*(1+jui-jli), & 
                 !!           (jli-2)*GRID_IHI_GC + i ) = 1.d0
                 ccAmat_2D( (jli-2)*GRID_IHI_GC + i, & 
                            (jli-2)*GRID_IHI_GC + i + GRID_IHI_GC*(1+jui-jli) ) = 1.d0
                 ccRHS_2D ( (jli-2)*GRID_IHI_GC + i ) = 0.d0
                 iNonZero = iNonZero + 2

                 !- kpd - Top cells left to right i
                 ccAmat_2D( (jui*GRID_IHI_GC) + i , &
                            (jui*GRID_IHI_GC) + i ) = -1.d0

                 !- kpd - Periodic BC
                 !!ccAmat_2D( (jui*GRID_IHI_GC) + i - GRID_IHI_GC*(1+jui-jli),   &
                 !!           (jui*GRID_IHI_GC) + i ) = 1.d0
                 ccAmat_2D( (jui*GRID_IHI_GC) + i, &
                            (jui*GRID_IHI_GC) + i - GRID_IHI_GC*(1+jui-jli) ) = 1.d0
                 ccRHS_2D ( (jui*GRID_IHI_GC) + i ) = 0.d0
                 iNonZero = iNonZero + 2

              end do

              !- kpd - Set Boundary Conditions for coeff matrix and RHS 
              do j=jli,jui
                 !- kpd - Left cells up and down j
                 ccAmat_2D( (j-1)*(GRID_IHI_GC)+(ili-1) , &
                            (j-1)*(GRID_IHI_GC)+(ili-1) ) = -1.d0                

                 !- kpd - Periodic BC
                 !!ccAmat_2D( (j-1)*(GRID_IHI_GC)+(ili-1) + (1+iui-ili), &
                 !!           (j-1)*(GRID_IHI_GC)+(ili-1) ) = 1.d0
                 ccAmat_2D( (j-1)*(GRID_IHI_GC)+(ili-1) , &
                            (j-1)*(GRID_IHI_GC)+(ili-1) + (1+iui-ili) ) = 1.d0
                 ccRHS_2D ( (j-1)*(GRID_IHI_GC)+(ili-1) ) = 0.d0
                 iNonZero = iNonZero + 2

                 !- kpd - Right cells up and down j
                 ccAmat_2D( (j-1)*(GRID_IHI_GC)+iui+1           , &
                            (j-1)*(GRID_IHI_GC)+iui+1 ) = -1.d0 

                 !- kpd - Periodic BC
                 !!ccAmat_2D( (j-1)*(GRID_IHI_GC)+iui+1 - (1+iui-ili) , &
                 !!           (j-1)*(GRID_IHI_GC)+iui+1 ) = 1.d0
                 ccAmat_2D( (j-1)*(GRID_IHI_GC)+iui+1 , &
                            (j-1)*(GRID_IHI_GC)+iui+1 - (1+iui-ili) ) = 1.d0
                 ccRHS_2D ( (j-1)*(GRID_IHI_GC)+iui+1 ) = 0.d0
                 iNonZero = iNonZero + 2

              end do 

              !- kpd - Another way to count Non-Zero Coeff Matrix values
              pIndex = 0
              do j = GRID_JLO_GC,GRID_IHI_GC*GRID_JHI_GC
                 do i = GRID_ILO_GC,GRID_IHI_GC*GRID_JHI_GC
                    if (ABS(ccAmat_2D(i,j)) .GT. 1e-8) then    
                       pIndex = pIndex + 1
                       !!- kpd - Sanity Check For Coeff Matrix
                       !print*,i,j,ccAmat_2D(i,j)
                    end if
                 end do
              end do
              print*,"pIndex",pIndex,iNonZero

              !- kpd - Put the coefficient matrix in compressed column format
              !call poisson_mg_ccFormat(iNonZero,GRID_ILO_GC,GRID_IHI_GC,ccAmat_2D,ccRHS_2D)
              call poisson_mg_ccFormat(pIndex,GRID_ILO_GC,GRID_IHI_GC,GRID_JLO_GC,GRID_JHI_GC,ccAmat_2D,ccRHS_2D,ccSol_2D)
            
              !- kpd - Set the final residual correction solution as "work"
              do j = jli, jui
                 do i = ili, iui
                    work(i,j,kli:kui,lb,1) = ccSol_2D(i,j)
                 end do
              end do

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unk,CENTER)
              !- kpd - Added for variable density...
              call Grid_releaseBlkPtr(lb,facexData,FACEX)
              call Grid_releaseBlkPtr(lb,faceyData,FACEY)

           endif

        enddo

     !************************************************************************
     !************************************************************************
     !************************************************************************

     !- kpd - 3-D Gauss Seidel
     else

        select case(gr_mgDiffOpDiscretize)

        case(2)    ! 2nd order Central Difference Solution
                   ! -------------------------------------

        do lb = 1, lnblocks
           if (lrefine(lb) == level) then
              call Grid_getBlkPhysicalSize(lb,size)

              delx = size(1)**2 * nxbinv
              dely = size(2)**2 * nybinv
              delz = size(3)**2 * nzbinv
              prefac = 0.5 / (dely*delz+delx*delz+delx*dely)


              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unk,CENTER)
              !- kpd - Added for variable density...
              call Grid_getBlkPtr(lb,facexData,FACEX)
              call Grid_getBlkPtr(lb,faceyData,FACEY)
              call Grid_getBlkPtr(lb,facezData,FACEZ)

              ksweep = 0
              do redblackpass = 1, 2
                 jsweep = ksweep
                 do k = kli, kui
                    isweep = jsweep
                    do j = jli, jui
                       do i = ili + isweep, iui, 2

                          !- kpd - Added for variable density...
                          MdensXL = ( facexData(RH1F_FACE_VAR,i  ,j  ,k  ) + &
                                      facexData(RH2F_FACE_VAR,i  ,j  ,k  ) )   ! Inverse density on left face.
                          MdensXR = ( facexData(RH1F_FACE_VAR,i+1,j  ,k  ) + &
                                      facexData(RH2F_FACE_VAR,i+1,j  ,k  ) )   ! Inverse density on right face.
                          MdensYL = ( faceyData(RH1F_FACE_VAR,i  ,j  ,k  ) + &
                                      faceyData(RH2F_FACE_VAR,i  ,j  ,k  ) )   ! Inverse density on bottom face.
                          MdensYR = ( faceyData(RH1F_FACE_VAR,i  ,j+1,k  ) + &
                                      faceyData(RH2F_FACE_VAR,i  ,j+1,k  ) )   ! Inverse density on top face.
                          MdensZL = ( facezData(RH1F_FACE_VAR,i  ,j  ,k  ) + &
                                      facezData(RH2F_FACE_VAR,i  ,j  ,k  ) )   ! Inverse density on front face.
                          MdensZR = ( facezData(RH1F_FACE_VAR,i  ,j  ,k+1) + &
                                      facezData(RH2F_FACE_VAR,i  ,j  ,k+1) )   ! Inverse density on back face.


                          !- kpd - Original Constant Density Implementation
                          !work(i,j,k,lb,1) = & 
                          !     &  prefac * & 
                          !     &  ( dely*delz*(work(i-1,j,k,lb,1) + work(i+1,j,k,lb,1)) + & 
                          !     &    delx*delz*(work(i,j-1,k,lb,1) + work(i,j+1,k,lb,1)) + & 
                          !     &    delx*dely*(work(i,j,k-1,lb,1) + work(i,j,k+1,lb,1)) - & 
                          !     &    delx*dely*delz*unk(irhs,i,j,k) )

                          !- kpd - Added for variable density...
                          !        div (1/rho) * (grad(P)) = dx^2 * (div(uStar + sig*K))
                          !        Note: unk(irhs,...) is the residual (r=RHS-A*P)
                          !        Note: work is the error correction (Ae=r, Pnew = Pold + e)
                          !        Note: delx is dx^2
                          work(i,j,k,lb,1) = &
                             ( MdensXR/delx * work(i+1,j,k,lb,1) +   &
                               MdensXL/delx * work(i-1,j,k,lb,1) +   &
                               MdensYR/dely * work(i,j+1,k,lb,1) +   &
                               MdensYL/dely * work(i,j-1,k,lb,1) +   &
                               MdensZR/delz * work(i,j,k+1,lb,1) +   &
                               MdensZL/delz * work(i,j,k-1,lb,1) -   &
                               unk(irhs,i,j,k)                   ) / &
                             ( MdensXR/delx + MdensXL/delx       +   &
                               MdensYR/dely + MdensYL/dely       +   &
                               MdensZR/delz + MdensZL/delz       )


                       enddo
                       isweep = 1 - isweep
                    enddo
                    jsweep = 1 - jsweep
                 enddo
                 ksweep = 1 - ksweep
              enddo
              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unk,CENTER)
              !- kpd - Added for variable density...
              call Grid_releaseBlkPtr(lb,facexData,FACEX)
              call Grid_releaseBlkPtr(lb,faceyData,FACEY)
              call Grid_releaseBlkPtr(lb,facezData,FACEZ)
           endif
        enddo

        case(4)    ! 4th order Central Difference Solution
                   ! -------------------------------------
        do lb = 1, lnblocks
           if (lrefine(lb) == level) then
              call Grid_getBlkPhysicalSize(lb,size)

              idelx = (size(1)**2 * nxbinv)**(-1.) !1/dx**2
              idely = (size(2)**2 * nybinv)**(-1.) !1/dy**2
              idelz = (size(3)**2 * nzbinv)**(-1.) !1/dz**2

              prefac = 576./1460.*(idelx+idely+idelz)**(-1.)    

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unk,CENTER)
              ksweep = 0
              do redblackpass = 1, 2
                 jsweep = ksweep
                 do k = kli, kui
                    isweep = jsweep
                    do j = jli, jui
                       do i = ili + isweep, iui, 2

                          work(i,j,k,lb,1) = & 
                               &  prefac * & 
                               &  (0.001736111111111 * &
                               &  (idelx*(     (work(i-3,j,k,lb,1) + work(i+3,j,k,lb,1)) - &
                               &           54.*(work(i-2,j,k,lb,1) + work(i+2,j,k,lb,1)) + &
                               &          783.*(work(i-1,j,k,lb,1) + work(i+1,j,k,lb,1)))+ &
                               &   idely*(     (work(i,j-3,k,lb,1) + work(i,j+3,k,lb,1)) - &
                               &           54.*(work(i,j-2,k,lb,1) + work(i,j+2,k,lb,1)) + &
                               &          783.*(work(i,j-1,k,lb,1) + work(i,j+1,k,lb,1)))+ &
                               &   idelz*(     (work(i,j,k-3,lb,1) + work(i,j,k+3,lb,1)) - &
                               &           54.*(work(i,j,k-2,lb,1) + work(i,j,k+2,lb,1)) + &
                               &          783.*(work(i,j,k-1,lb,1) + work(i,j,k+1,lb,1)))) &
                               &   - unk(irhs,i,j,k) )

                       enddo
                       isweep = 1 - isweep
                    enddo
                    jsweep = 1 - jsweep
                 enddo
                 ksweep = 1 - ksweep
              enddo
              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unk,CENTER)
           endif
        enddo

        end select


     endif

     ! Now compute the error norm (change in LHS from one iteration
     ! to the next).  This will be used to limit the number of
     ! iterations (see above).  Here we use the maximum norm rather
     ! than the L2 norm computed by mg_norm.

     ! However, don't start checking the convergence criterion
     ! until we've done ten iterations.  This saves work when doing
     ! the regular smoothing on fine levels, which isn't expected
     ! to produce convergence by itself.


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

       !**************************************************************
       !print*,"Iterative error norms: ",lerrnorm(1),lerrnorm(2),error
       !**************************************************************

     endif 

     !End of loop.  Iterate until nsmooth iterations have been
     !    performed, or until the convergence criterion has been met.
     done = (iter == nsmooth) .or. ((iter > iterating_to_convergence_limit) .and. (error > 0.))

     !call timer_stop("relax calc")

     ! fill work's guardcells
     if (.not. done) &
     call gr_mgBndry (level, ilhs, nguard_work, 0, MG_EXCHANGE_WORK, MG_CONTINUE_SERIES) 

  enddo

!================================================================================================
!================================================================================================

  !***********************************************************
  !print*,"Final Iterative Error",iter,lerrnorm(1),lerrnorm(2)
  !***********************************************************

  ! *** Why is this needed ?  It is called above !!!
  call gr_mgBndry(level, ilhs, nguard_work, 0, MG_EXCHANGE_WORK, MG_END_SERIES)

  !-----------------------------------------------------------------------
  !- kpd - Store the final solution (work) in ilhs
  !        Note: work is only the error correction (Ae=r, Pnew = Pold + e)
  !-----------------------------------------------------------------------
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

  !===============================================================================

  !call timer_stop("relax")

  return

end subroutine poisson_mg_relax
