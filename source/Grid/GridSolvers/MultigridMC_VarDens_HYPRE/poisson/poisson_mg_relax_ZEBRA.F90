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


subroutine poisson_mg_relax_ZEBRA (level, irhs, ilhs, nsmooth, idenvar, levelmax)

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

  integer       :: level, irhs, ilhs, nsmooth, idebug, idenvar, levelmax
  integer       :: firstpass = 1

  integer       :: iter, i, j, k, lb, ierr, Neff, lnblocks,iii
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
                                       facevarx,facevary,facevarz
  !- kpd - 
  real :: MdensXL, MdensXR, MdensYL, MdensYR, MdensZL, MdensZR

  real, dimension(NXB) :: ax,bx,cx,fx
  real, dimension(NYB) :: ay,by,cy,fy
  real, dimension(NZB) :: az,bz,cz,fz

  integer, parameter :: nxb = NXB
  integer, parameter :: nyb = NYB
  integer, parameter :: nzb = NZB
  integer, parameter :: ndim = NDIM


  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox

  integer blockcount,ii,jj,blockID,i1,j1,k1,zebrapass

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

  !- KPD - !!!!!!!!!!!!!!!!!!!!
  !- KPD - !!!!!!!!!!!!!!!!!!!!
  !- KPD - !!!!!!!!!!!!!!!!!!!!
  !- kpd - I put nodetype restoration in b/c paramesh was crashing for AMR with >1 level
  call mg_restore_nodetypes (level)
  call gr_mgBndry (level, ilhs, 2, 0, MG_COPY_UNK_TO_WORK, MG_BEGIN_SERIES) 

  iii=0

  do while ( .not. done )

     iii = iii + 1

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
              call Grid_getBlkPtr(lb,facevarx ,FACEX)
              call Grid_getBlkPtr(lb,facevary ,FACEY)

              !if (iii .eq. 1) then
              !print*,"poisson_mg_relax_RBGS.F90 Blocks:",myPE,lb
              !end if

              !- kpd - These are all 1/dx^2 or 1/dy^2
              delx = size(1)**2 * nxbinv
              dely = size(2)**2 * nybinv
              !prefac = 0.5 / (delx+dely)

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unk,CENTER)


              ! Do odd-even pass on x-lines              
              do zebrapass = 0, 1
                 do k = kli,kui
                   do j = jli+zebrapass, jui, 2
                      i1 = 0
                      do i = ili,iui

                       !---------------------------------------------------------------
                       !- kpd - Added for variable density. Calculate face densities...
                       !---------------------------------------------------------------
!                       if (level == levelmax) then
                          MdensXL = ( facexData(RH1F_FACE_VAR,i  ,j  ,k) + &
                                      facexData(RH2F_FACE_VAR,i  ,j  ,k) )   ! Inverse density on left face.
                          MdensXR = ( facexData(RH1F_FACE_VAR,i+1,j  ,k) + &
                                      facexData(RH2F_FACE_VAR,i+1,j  ,k) )   ! Inverse density on right face.
                          MdensYL = ( faceyData(RH1F_FACE_VAR,i  ,j  ,k) + &
                                      faceyData(RH2F_FACE_VAR,i  ,j  ,k) )   ! Inverse density on bottom face.
                          MdensYR = ( faceyData(RH1F_FACE_VAR,i  ,j+1,k) + &
                                      faceyData(RH2F_FACE_VAR,i  ,j+1,k) )   ! Inverse density on top face.
!                       else
!                          MdensXL = facevarx(idenvar,i  ,j  ,k)
!                          MdensXR = facevarx(idenvar,i+1,j  ,k)
!                          MdensYL = facevary(idenvar,i  ,j  ,k)
!                          MdensYR = facevary(idenvar,i  ,j+1,k)
!                       end if

                       i1 = i1 + 1

                       ! Coefficients of a, b, c:
                       ax(i1) =   MdensXR/delx + MdensXL/delx + MdensYR/dely + MdensYL/dely
                       bx(i1) = - MdensXL/delx
                       cx(i1) = - MdensXR/delx

                       ! Right hand side:
                       fx(i1) = MdensYR/dely * work(i,j+1,k,lb,1) +   &
                                MdensYL/dely * work(i,j-1,k,lb,1) -   &
                                unk(irhs,i,j,k) 

                       ! Boundary conditions:
                       if (i1 == 1) then
                          bx(i1) = 0.
                          fx(i1) = fx(i1) + MdensXL/delx * work(i-1,j,k,lb,1)
                       elseif (i1 == NXB) then 
                          cx(i1) = 0.
                          fx(i1) = fx(i1) + MdensXR/delx * work(i+1,j,k,lb,1)
                       end if

                      enddo

                      ! Call tridiagonal solver:
                      call solve_tridiag(bx,ax,cx,fx,work(ili:iui,j,k,lb,1),NXB)

                   enddo
                enddo
             enddo

             ! Do odd-even pass on y-lines              
             do zebrapass = 0, 1
                do k=kli,kui
                   do i = ili+zebrapass, iui, 2
                      j1 = 0
                      do j = jli,jui

                       !---------------------------------------------------------------
                       !- kpd - Added for variable density. Calculate face densities...
                       !---------------------------------------------------------------
!                       if (level == levelmax) then
                          MdensXL = ( facexData(RH1F_FACE_VAR,i  ,j  ,k) + &
                                      facexData(RH2F_FACE_VAR,i  ,j  ,k) )   ! Inverse density on left face.
                          MdensXR = ( facexData(RH1F_FACE_VAR,i+1,j  ,k) + &
                                      facexData(RH2F_FACE_VAR,i+1,j  ,k) )   ! Inverse density on right face.
                          MdensYL = ( faceyData(RH1F_FACE_VAR,i  ,j  ,k) + &
                                      faceyData(RH2F_FACE_VAR,i  ,j  ,k) )   ! Inverse density on bottom face.
                          MdensYR = ( faceyData(RH1F_FACE_VAR,i  ,j+1,k) + &
                                      faceyData(RH2F_FACE_VAR,i  ,j+1,k) )   ! Inverse density on top face.
!                       else
!                          MdensXL = facevarx(idenvar,i  ,j  ,k)
!                          MdensXR = facevarx(idenvar,i+1,j  ,k)
!                          MdensYL = facevary(idenvar,i  ,j  ,k)
!                          MdensYR = facevary(idenvar,i  ,j+1,k)
!                       end if

                       j1 = j1 + 1

                       ! Coefficients of a, b, c:
                       ay(j1) =   MdensXR/delx + MdensXL/delx + MdensYR/dely + MdensYL/dely
                       by(j1) = - MdensYL/dely
                       cy(j1) = - MdensYR/dely

                       ! Right hand side:
                       fy(j1) =  MdensXR/delx * work(i+1,j,k,lb,1) +   &
                                 MdensXL/delx * work(i-1,j,k,lb,1) -   &
                                 unk(irhs,i,j,k) 

                       ! Boundary conditions:
                       if (j1 == 1) then
                          by(j1) = 0.
                          fy(j1) = fy(j1) + MdensYL/dely * work(i,j-1,k,lb,1) 
                       elseif (j1 == NYB) then 
                          cy(j1) = 0.
                          fy(j1) = fy(j1) + MdensYR/dely * work(i,j+1,k,lb,1)
                       end if

                      enddo

                      ! Call tridiagonal solver:
                      call solve_tridiag(by,ay,cy,fy,work(i,jli:jui,k,lb,1),NYB)

                   enddo
                enddo
              enddo



              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unk,CENTER)
              call Grid_releaseBlkPtr(lb,facexData,FACEX)
              call Grid_releaseBlkPtr(lb,faceyData,FACEY)
              call Grid_releaseBlkPtr(lb,facevarx ,FACEX)
              call Grid_releaseBlkPtr(lb,facevary ,FACEY)

           endif
        enddo


     else    !ndim == 3

        select case(gr_mgDiffOpDiscretize)

        case(2)                ! 2nd order Central Difference Solution
        do lb = 1, lnblocks    !--------------------------------------
           if (lrefine(lb) == level) then
              call Grid_getBlkPhysicalSize(lb,size)

              delx = size(1)**2 * nxbinv
              dely = size(2)**2 * nybinv
              delz = size(3)**2 * nzbinv
              prefac = 0.5 / (dely*delz+delx*delz+delx*dely)

              !- kpd
              call Grid_getBlkPtr(lb,facexData,FACEX)
              call Grid_getBlkPtr(lb,faceyData,FACEY)
              call Grid_getBlkPtr(lb,facezData,FACEZ)
              call Grid_getBlkPtr(lb,facevarx ,FACEX)
              call Grid_getBlkPtr(lb,facevary ,FACEY)
              call Grid_getBlkPtr(lb,facevarz ,FACEZ)

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unk,CENTER)


              ! Do odd-even pass on x-lines              
              do zebrapass = 0, 1
                 do k = kli,kui
                   do j = jli+zebrapass, jui, 2
                      i1 = 0
                      do i = ili,iui


                       !---------------------------------------------------------------
                       !- kpd - Added for variable density. Calculate face densities...
                       !---------------------------------------------------------------
!                       if (level == levelmax) then
                          MdensXL = ( facexData(RH1F_FACE_VAR,i  ,j  ,k) + &
                                      facexData(RH2F_FACE_VAR,i  ,j  ,k) )   ! Inverse density on left face.
                          MdensXR = ( facexData(RH1F_FACE_VAR,i+1,j  ,k) + &
                                      facexData(RH2F_FACE_VAR,i+1,j  ,k) )   ! Inverse density on right face.
                          MdensYL = ( faceyData(RH1F_FACE_VAR,i  ,j  ,k) + &
                                      faceyData(RH2F_FACE_VAR,i  ,j  ,k) )   ! Inverse density on bottom face.
                          MdensYR = ( faceyData(RH1F_FACE_VAR,i  ,j+1,k) + &
                                      faceyData(RH2F_FACE_VAR,i  ,j+1,k) )   ! Inverse density on top face.
                          MdensZL = ( facezData(RH1F_FACE_VAR,i  ,j  ,k) + &
                                      facezData(RH2F_FACE_VAR,i  ,j  ,k) )   ! Inverse density on back face.
                          MdensZR = ( facezData(RH1F_FACE_VAR,i  ,j,k+1) + &
                                      facezData(RH2F_FACE_VAR,i  ,j,k+1) )   ! Inverse density on front face.
!                       else
!                          MdensXL = facevarx(idenvar,i  ,j  ,k)
!                          MdensXR = facevarx(idenvar,i+1,j  ,k)
!                          MdensYL = facevary(idenvar,i  ,j  ,k)
!                          MdensYR = facevary(idenvar,i  ,j+1,k)
!                          MdensZL = facevarz(idenvar,i  ,j  ,k)
!                          MdensZR = facevarz(idenvar,i  ,j,k+1)
!                       end if

                       i1 = i1 + 1

                       ! Coefficients of a, b, c:
                       ax(i1) =   MdensXR/delx + MdensXL/delx + &
                                  MdensYR/dely + MdensYL/dely + &
                                  MdensZR/delz + MdensZL/delz 
                       bx(i1) = - MdensXL/delx
                       cx(i1) = - MdensXR/delx

                       ! Right hand side:
                       fx(i1) = MdensYR/dely * work(i,j+1,k,lb,1) +   &
                                MdensYL/dely * work(i,j-1,k,lb,1) +   &
                                MdensZR/delz * work(i,j,k+1,lb,1) +   &
                                MdensZL/delz * work(i,j,k-1,lb,1) -   &
                                unk(irhs,i,j,k) 

                       ! Boundary conditions:
                       if (i1 == 1) then
                          bx(i1) = 0.
                          fx(i1) = fx(i1) + MdensXL/delx * work(i-1,j,k,lb,1)
                       elseif (i1 == NXB) then 
                          cx(i1) = 0.
                          fx(i1) = fx(i1) + MdensXR/delx * work(i+1,j,k,lb,1)
                       end if

                      enddo

                      ! Call tridiagonal solver:
                      call solve_tridiag(bx,ax,cx,fx,work(ili:iui,j,k,lb,1),NXB)

                   enddo
                enddo
             enddo



             ! Do odd-even pass on y-lines              
             do zebrapass = 0, 1
                do k=kli+zebrapass, kui, 2
                   do i = ili,iui
                      j1 = 0
                      do j = jli,jui

                       !---------------------------------------------------------------
                       !- kpd - Added for variable density. Calculate face densities...
                       !---------------------------------------------------------------
!                       if (level == levelmax) then
                          MdensXL = ( facexData(RH1F_FACE_VAR,i  ,j  ,k) + &
                                      facexData(RH2F_FACE_VAR,i  ,j  ,k) )   ! Inverse density on left face.
                          MdensXR = ( facexData(RH1F_FACE_VAR,i+1,j  ,k) + &
                                      facexData(RH2F_FACE_VAR,i+1,j  ,k) )   ! Inverse density on right face.
                          MdensYL = ( faceyData(RH1F_FACE_VAR,i  ,j  ,k) + &
                                      faceyData(RH2F_FACE_VAR,i  ,j  ,k) )   ! Inverse density on bottom face.
                          MdensYR = ( faceyData(RH1F_FACE_VAR,i  ,j+1,k) + &
                                      faceyData(RH2F_FACE_VAR,i  ,j+1,k) )   ! Inverse density on top face.
                          MdensZL = ( facezData(RH1F_FACE_VAR,i  ,j  ,k) + &
                                      facezData(RH2F_FACE_VAR,i  ,j  ,k) )   ! Inverse density on back face.
                          MdensZR = ( facezData(RH1F_FACE_VAR,i  ,j,k+1) + &
                                      facezData(RH2F_FACE_VAR,i  ,j,k+1) )   ! Inverse density on front face.
!                       else
!                          MdensXL = facevarx(idenvar,i  ,j  ,k)
!                          MdensXR = facevarx(idenvar,i+1,j  ,k)
!                          MdensYL = facevary(idenvar,i  ,j  ,k)
!                          MdensYR = facevary(idenvar,i  ,j+1,k)
!                          MdensZL = facevarz(idenvar,i  ,j  ,k)
!                          MdensZR = facevarz(idenvar,i  ,j,k+1)
!                       end if

                       j1 = j1 + 1

                       ! Coefficients of a, b, c:
                       ay(j1) =   MdensXR/delx + MdensXL/delx + &
                                  MdensYR/dely + MdensYL/dely + &
                                  MdensZR/delz + MdensZL/delz 
                       by(j1) = - MdensYL/dely
                       cy(j1) = - MdensYR/dely

                       ! Right hand side:
                       fy(j1) =  MdensXR/delx * work(i+1,j,k,lb,1) +   &
                                 MdensXL/delx * work(i-1,j,k,lb,1) +   &
                                 MdensZR/delz * work(i,j,k+1,lb,1) +   &
                                 MdensZL/delz * work(i,j,k-1,lb,1) -   &
                                 unk(irhs,i,j,k) 

                       ! Boundary conditions:
                       if (j1 == 1) then
                          by(j1) = 0.
                          fy(j1) = fy(j1) + MdensYL/dely * work(i,j-1,k,lb,1) 
                       elseif (j1 == NYB) then 
                          cy(j1) = 0.
                          fy(j1) = fy(j1) + MdensYR/dely * work(i,j+1,k,lb,1)
                       end if

                      enddo

                      ! Call tridiagonal solver:
                      call solve_tridiag(by,ay,cy,fy,work(i,jli:jui,k,lb,1),NYB)

                   enddo
                enddo
              enddo


             ! Do odd-even pass on z-lines              
             do zebrapass = 0, 1
                do j = jli,jui 
                   do i = ili+zebrapass, iui, 2
                      k1 = 0
                      do k=kli,kui

                       !---------------------------------------------------------------
                       !- kpd - Added for variable density. Calculate face densities...
                       !---------------------------------------------------------------
!                       if (level == levelmax) then
                          MdensXL = ( facexData(RH1F_FACE_VAR,i  ,j  ,k) + &
                                      facexData(RH2F_FACE_VAR,i  ,j  ,k) )   ! Inverse density on left face.
                          MdensXR = ( facexData(RH1F_FACE_VAR,i+1,j  ,k) + &
                                      facexData(RH2F_FACE_VAR,i+1,j  ,k) )   ! Inverse density on right face.
                          MdensYL = ( faceyData(RH1F_FACE_VAR,i  ,j  ,k) + &
                                      faceyData(RH2F_FACE_VAR,i  ,j  ,k) )   ! Inverse density on bottom face.
                          MdensYR = ( faceyData(RH1F_FACE_VAR,i  ,j+1,k) + &
                                      faceyData(RH2F_FACE_VAR,i  ,j+1,k) )   ! Inverse density on top face.
                          MdensZL = ( facezData(RH1F_FACE_VAR,i  ,j  ,k) + &
                                      facezData(RH2F_FACE_VAR,i  ,j  ,k) )   ! Inverse density on back face.
                          MdensZR = ( facezData(RH1F_FACE_VAR,i  ,j,k+1) + &
                                      facezData(RH2F_FACE_VAR,i  ,j,k+1) )   ! Inverse density on front face.
!                       else
!                          MdensXL = facevarx(idenvar,i  ,j  ,k)
!                          MdensXR = facevarx(idenvar,i+1,j  ,k)
!                          MdensYL = facevary(idenvar,i  ,j  ,k)
!                          MdensYR = facevary(idenvar,i  ,j+1,k)
!                          MdensZL = facevarz(idenvar,i  ,j  ,k)
!                          MdensZR = facevarz(idenvar,i  ,j,k+1)
!                       end if

                       k1 = k1 + 1

                       ! Coefficients of a, b, c:
                       az(k1) =   MdensXR/delx + MdensXL/delx + &
                                  MdensYR/dely + MdensYL/dely + &
                                  MdensZR/delz + MdensZL/delz 
                       bz(k1) = - MdensZL/delz
                       cz(k1) = - MdensZR/delz

                       ! Right hand side:
                       fz(k1) =  MdensXR/delx * work(i+1,j,k,lb,1) +   &
                                 MdensXL/delx * work(i-1,j,k,lb,1) +   &
                                 MdensYR/dely * work(i,j+1,k,lb,1) +   &
                                 MdensYL/dely * work(i,j-1,k,lb,1) -   &
                                 unk(irhs,i,j,k) 

                       ! Boundary conditions:
                       if (k1 == 1) then
                          bz(k1) = 0.
                          fz(k1) = fz(k1) + MdensZL/delz * work(i,j,k-1,lb,1)
                       elseif (k1 == NZB) then 
                          cz(k1) = 0.
                          fz(k1) = fz(k1) + MdensZR/delz * work(i,j,k+1,lb,1)
                       end if

                      enddo

                      ! Call tridiagonal solver:
                      call solve_tridiag(bz,az,cz,fz,work(i,j,kli:kui,lb,1),NZB)

                   enddo
                enddo
              enddo

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unk,CENTER)
              call Grid_releaseBlkPtr(lb,facexData,FACEX)
              call Grid_releaseBlkPtr(lb,faceyData,FACEY)
              call Grid_releaseBlkPtr(lb,facezData,FACEZ)
              call Grid_releaseBlkPtr(lb,facevarx ,FACEX)
              call Grid_releaseBlkPtr(lb,facevary ,FACEY)
              call Grid_releaseBlkPtr(lb,facevarz ,FACEZ)
           endif
        enddo

        case(4)    ! 4th order Central difference Solution

           print*,"ERROR: 4th Order Central Not Implemented For Variable Density"

        end select

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
! Q -  Why is this needed ?  It is called above !!! kpd
! A -  Because BC for last iteration is not satisfied in the loop !!! Akash

!!!  call mg_bndry(level, ilhs, nguard_work, 0, EXCHANGE_WORK, CONTINUE_SERIES)
  call gr_mgBndry(level, ilhs, nguard_work, 0, MG_EXCHANGE_WORK, MG_END_SERIES)

  do lb = 1, lnblocks
     if (lrefine(lb) == level) then

        ! Point to blocks center vars:
        call Grid_getBlkPtr(lb,unk,CENTER)

        !print*,"RELrbgs LB",lb

        do k = kli, kui
           do j = jli, jui
              do i = ili, iui
                 unk(ilhs,i,j,k) = work(i,j,k,lb,1)

!if (lb .eq. 58 .or. lb .eq. 57 ) then
!   print*,lb,i,j,unk(ilhs,i,j,k)
!end if

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
end subroutine poisson_mg_relax_ZEBRA


!---------------------------------------------------------------------------------

     ! Tridiagonal solver Routine from wiki:
     subroutine solve_tridiag(a,b,c,d,x,n)
      implicit none
!        a - sub-diagonal (means it is the diagonal below the main diagonal)
!        b - the main diagonal
!        c - sup-diagonal (means it is the diagonal above the main diagonal)
!        d - right part
!        x - the answer
!        n - number of equations
 
        integer,intent(in) :: n
        real,dimension(n),intent(in) :: a,b,c,d
        real,dimension(n),intent(out) :: x
        real,dimension(n) :: cp,dp
        real :: m
        integer i
 
! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
         do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/m
           dp(i) = (d(i)-dp(i-1)*a(i))/m
         enddo
! initialize x
         x(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)
        end do
 
    end subroutine solve_tridiag
