!*******************************************************************************

!  Routine:     mg_residual()

!  Description: Compute the residual of the equation to be solved using
!               multigrid, given a guess at the solution on a particular
!               level and the source term on that level.  This version
!               implements residuals for the Poisson equation.

!  Parameters:  level       Level to compute the residual on.
!               irhs        Right-hand side (source) of equation.
!               ilhs        Left-hand side (solution) of equation.
!               ires        Index of variable to receive the residual.
!               leaf_only   If nonzero, only compute residual on leaf nodes
!                             on the specified level.
!               flux_cons   If nonzero, match fluxes at boundaries with
!                             finer blocks


subroutine poisson_mg_residual (level, irhs, ilhs, ires, leaf_only, & 
     &                        flux_cons)

  !===============================================================================

  use mg_common

  use Grid_data, ONLY : gr_meshMe,gr_meshNumProcs

  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getBlkPtr,       &
                                Grid_releaseBlkPtr,   &
                                Grid_getBlkPhysicalSize, &
                                Grid_getListOfBlocks, &
                                Grid_getBlkIndexLimits


  use physicaldata, only : flux_x,flux_y,flux_z
  use workspace, ONLY: work
  use tree, only : nodetype,lrefine
  use paramesh_dimensions, only : nguard_work
  use paramesh_interfaces, only :  amr_flux_conserve


  implicit none
#include "constants.h"
#include "Multigrid.h"
#include "Flash.h"  
#include "mpif.h"

  integer :: level, irhs, ilhs, ires, leaf_only, flux_cons

  integer :: i, j, k, lb, ierr
  real    :: delx, dely, delz
  real    :: norm
  logical :: some_proc_needs_updating, my_proc_needs_updating
  logical, dimension(MAXBLOCKS) :: update_fluxes
  logical, dimension(MAXBLOCKS) :: update_values

  integer, parameter       :: MAXDIMS2 = 3
  integer, parameter       :: MFACES2 = 2*MAXDIMS2
  real, dimension(MAXDIMS2) :: size
  real, dimension(MFACES2)  :: neigh2, neigh_type



  real    :: mgfluxx(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC, & 
       &                   MAXBLOCKS)
  real    :: mgfluxy(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC, & 
       &                   MAXBLOCKS)
  real    :: mgfluxz(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC, & 
       &                   MAXBLOCKS)

  real, pointer, save, dimension(:,:,:,:,:) :: flux_x_ptr, flux_y_ptr, flux_z_ptr

  real, pointer, save, dimension(:,:,:,:) :: unkt

  !- kpd -
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData
  !- kpd - 
  real :: MdensXL, MdensXR, MdensYL, MdensYR, MdensZL, MdensZR

  integer :: lnblocks2
  integer, save :: myPE2,num_PEs

  real dx,dy,dz,dxdz,dydz,dxdy

  integer, parameter :: nxb = NXB
  integer, parameter :: nyb = NYB
  integer, parameter :: nzb = NZB
  integer, parameter :: nguard = NGUARD
  integer, parameter :: ndim = NDIM

  logical, save :: first_call = .true.

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox

  integer blockCount,ii,jj,blockID

  real bsize(MDIM),coord(MDIM)

  integer blockList(MAXBLOCKS)

  !=============================================================================

!  call timer_start("residual")

!  if (leaf_only.ne.0) then
!     call timer_start("leaf_only")
!  else
!     call timer_start("not leaf_only")
!  end if

  if (first_call) then
     first_call = .false.
     myPE2 = gr_meshMe
     num_PEs = gr_meshNumProcs
     flux_x_ptr => flux_x
     flux_y_ptr => flux_y
     flux_z_ptr => flux_z
  end if
  call Grid_getLocalNumBlks(lnblocks2)


  flux_x_ptr(FLXMC_FLUX,:,:,:,:) = 0.
  flux_y_ptr(FLXMC_FLUX,:,:,:,:) = 0.
  flux_z_ptr(FLXMC_FLUX,:,:,:,:) = 0.

  my_proc_needs_updating = .false.
  some_proc_needs_updating = .false.

  do lb = 1, lnblocks2
     update_values(lb) = (lrefine(lb) == level)
     if (leaf_only /= 0) then
        update_values(lb) = update_values(lb) .and. (nodetype_save(lb) == 1)
     end if
     if (update_values(lb)) then
        my_proc_needs_updating = .true.
     end if
  end do

  call mpi_allreduce(my_proc_needs_updating, some_proc_needs_updating, &
       1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)

  if (some_proc_needs_updating) then
     do lb = 1, lnblocks2
        update_fluxes(lb) = (lrefine(lb) == level)
        if (flux_cons /= 0) then 
           update_fluxes(lb) = update_fluxes(lb) .or. (lrefine(lb) == level+1)
        end if
        if (leaf_only /= 0) then
           update_fluxes(lb) = update_fluxes(lb) .and. (nodetype_save(lb) == 1)
        end if
        if (update_fluxes(lb)) then
        end if
     end do
     

     !               Update boundary zones of the LHS array.
     !               Currently the mesh package doesn't supply a general boundary
     !               update routine, so we must copy the LHS into the "work" array,
     !               then update its boundaries.  This is handled by mg_bndry().

! ***
     call gr_mgBndry (level, ilhs, nguard_work, leaf_only, MG_COPY_UNK_TO_WORK, MG_STANDALONE) !MG_CONTINUE_SERIES

     !               Compute x, y, and z "fluxes."  Copy block boundary fluxes to
     !               the arrays used by the PARAMESH flux conservation routines.
     if (ndim == 1) then

        do lb = 1, lnblocks2

           if (update_fluxes(lb)) then
              

              ! Get BlockSize:
              call Grid_getBlkPhysicalSize(lb,size)
              delx = nxb/size(1)
              do i = ili-1, iui
                 mgfluxx(i,jli:jui,kli:kui,lb) = & 
                      &          delx * (work(i+1,jli:jui,kli:kui,lb,1) - & 
                      &                  work(i,jli:jui,kli:kui,lb,1))
              enddo
           endif

           flux_x_ptr(FLXMC_FLUX,1,:,:,lb) = mgfluxx(ili-1,jli:jui,kli:kui,lb)
           flux_x_ptr(FLXMC_FLUX,2,:,:,lb) = mgfluxx(iui,jli:jui,kli:kui,lb)



        enddo

     elseif (ndim == 2) then

        do lb = 1, lnblocks2

           ! Get BlockSize:
           call Grid_getBlkPhysicalSize(lb,size)
           !- kpd - 
           call Grid_getBlkPtr(lb,facexData,FACEX)
           call Grid_getBlkPtr(lb,faceyData,FACEY)

           dx = size(1)/real(nxb)
           dy = size(2)/real(nyb)


           if (update_fluxes(lb)) then

              !print*,"Block List inside poisson_mg_residual.F90. lb=",lb,"blockCount=",lnblocks2

              delx = real(nxb)/size(1)
              dely = real(nyb)/size(2)
              do j = jli-1, jui
                 do i = ili-1, iui

                    !- kpd - Added for variable density. Calculate face densities...
                    MdensXR = ( facexData(RH1F_FACE_VAR,i+1,j  ,1) + &
                                facexData(RH2F_FACE_VAR,i+1,j  ,1) )   ! Inverse density on right face.
                    MdensYR = ( faceyData(RH1F_FACE_VAR,i  ,j+1,1) + &
                                faceyData(RH2F_FACE_VAR,i  ,j+1,1) )   ! Inverse density on top face.`

                    if (MdensXR .le. 0. .OR. MdensYR .le. 0. ) then
                       print*,"ERROR 1 poisson_mg_residual.F90: Density Equals Zero at a FACE",lb,i,j,MdensXR,MdensYR
                    elseif (MdensXR.gt.1000.0 .OR. MdensYR.gt.1000.0 ) then
                       print*,"ERROR 1 poisson_mg_residual.F90: Inverse Density Greater than 1.0",lb,i,j,MdensXR,MdensYR
                    end if

                    !---------------------------------------
                    !- kpd - Constant density implementation
                    !---------------------------------------
                    !mgfluxx(i,j,kli:kui,lb) = & 
                    !     &            delx * (work(i+1,j,kli:kui,lb,1) - & 
                    !     &                    work(i,j,kli:kui,lb,1))
                    !mgfluxy(i,j,kli:kui,lb) = & 
                    !     &            dely * (work(i,j+1,kli:kui,lb,1) - & 
                    !     &                    work(i,j,kli:kui,lb,1))

                    !- kpd - Variable density implementation
                    !---------------------------------------
                    mgfluxx(i,j,kli:kui,lb) = & 
                                MdensXR * delx * (work(i+1,j,kli:kui,lb,1) - & 
                                                  work(i  ,j,kli:kui,lb,1))
                    mgfluxy(i,j,kli:kui,lb) = & 
                                MdensYR * dely * (work(i,j+1,kli:kui,lb,1) - & 
                                                  work(i,j  ,kli:kui,lb,1))

                 enddo
              enddo
           endif

           flux_x_ptr(FLXMC_FLUX,1,:,:,lb) = mgfluxx(ili-1,jli:jui,kli:kui,lb)*dy
           flux_x_ptr(FLXMC_FLUX,2,:,:,lb) = mgfluxx(iui,jli:jui,kli:kui,lb)*dy
           flux_y_ptr(FLXMC_FLUX,:,1,:,lb) = mgfluxy(ili:iui,jli-1,kli:kui,lb)*dx
           flux_y_ptr(FLXMC_FLUX,:,2,:,lb) = mgfluxy(ili:iui,jui,kli:kui,lb)*dx

           !- kpd - 
           call Grid_releaseBlkPtr(lb,facexData,FACEX)
           call Grid_releaseBlkPtr(lb,faceyData,FACEY)

        enddo

     else

        select case(gr_mgDiffOpDiscretize)

        case(2)  ! 2nd order Central difference Solution
        do lb = 1, lnblocks2

           ! Get BlockSize:
           call Grid_getBlkPhysicalSize(lb,size)
           !- kpd - 
           call Grid_getBlkPtr(lb,facexData,FACEX)
           call Grid_getBlkPtr(lb,faceyData,FACEY)
           call Grid_getBlkPtr(lb,facezData,FACEZ)

           dxdy =(size(1)*real(nxb)**(-1.)) * (size(2)*real(nyb)**(-1.))
           dydz =(size(2)*real(nyb)**(-1.)) * (size(3)*real(nzb)**(-1.))
           dxdz =(size(1)*real(nxb)**(-1.)) * (size(3)*real(nzb)**(-1.)) 

           if (update_fluxes(lb)) then

              delx = nxb/size(1)
              dely = nyb/size(2)
              delz = nzb/size(3)
              do k = kli-1, kui
                 do j = jli-1, jui
                    do i = ili-1, iui

                       !- kpd - Added for variable density. Calculate face densities...
                       MdensXR = ( facexData(RH1F_FACE_VAR,i+1,j  ,k) + &
                                   facexData(RH2F_FACE_VAR,i+1,j  ,k) )   ! Inverse density on right face.
                       MdensYR = ( faceyData(RH1F_FACE_VAR,i  ,j+1,k) + &
                                   faceyData(RH2F_FACE_VAR,i  ,j+1,k) )   ! Inverse density on top face.
                       MdensZR = ( facezData(RH1F_FACE_VAR,i  ,j,k+1) + &
                                   facezData(RH2F_FACE_VAR,i  ,j,k+1) )   ! Inverse density on front face.

                       if (MdensXR.le.0. .OR. MdensYR.le.0. .OR. MdensZR.le.0.) then
                          print*,"ERROR 2 poisson_mg_residual.F90: Density Equals Zero at a FACE",lb,i,j,MdensXR,MdensYR
                       elseif (MdensXR.gt.1000.0 .OR. MdensYR.gt.1000.0 .OR. MdensZR.gt.1000.0) then
                          print*,"ERROR 2 poisson_mg_residual.F90: Inverse Density Greater than 1.0",lb,i,j,MdensXR,MdensYR
                       end if

                       !---------------------------------------
                       !- kpd - Variable density implementation
                       !---------------------------------------
                       mgfluxx(i,j,k,lb) = & 
                            &              MdensXR * delx * (work(i+1,j,k,lb,1) - work(i,j,k,lb,1))
                       mgfluxy(i,j,k,lb) = & 
                            &              MdensYR * dely * (work(i,j+1,k,lb,1) - work(i,j,k,lb,1))
                       mgfluxz(i,j,k,lb) = & 
                            &              MdensZR * delz * (work(i,j,k+1,lb,1) - work(i,j,k,lb,1))

                    enddo
                 enddo
              enddo
           endif
           flux_x_ptr(FLXMC_FLUX,1,:,:,lb) = mgfluxx(ili-1,jli:jui,kli:kui,lb)*dydz
           flux_x_ptr(FLXMC_FLUX,2,:,:,lb) = mgfluxx(iui,jli:jui,kli:kui,lb)*dydz
           flux_y_ptr(FLXMC_FLUX,:,1,:,lb) = mgfluxy(ili:iui,jli-1,kli:kui,lb)*dxdz
           flux_y_ptr(FLXMC_FLUX,:,2,:,lb) = mgfluxy(ili:iui,jui,kli:kui,lb)*dxdz
           flux_z_ptr(FLXMC_FLUX,:,:,1,lb) = mgfluxz(ili:iui,jli:jui,kli-1,lb)*dxdy
           flux_z_ptr(FLXMC_FLUX,:,:,2,lb) = mgfluxz(ili:iui,jli:jui,kui,lb)*dxdy

           !- kpd - 
           call Grid_releaseBlkPtr(lb,facexData,FACEX)
           call Grid_releaseBlkPtr(lb,faceyData,FACEY)
           call Grid_releaseBlkPtr(lb,facezData,FACEZ)

        enddo

        case(4)    ! 4th order Central difference Solution
        do lb = 1, lnblocks2

           ! Get BlockSize:
           call Grid_getBlkPhysicalSize(lb,size)

           dxdy =(size(1)*real(nxb)**(-1.)) * (size(2)*real(nyb)**(-1.))
           dydz =(size(2)*real(nyb)**(-1.)) * (size(3)*real(nzb)**(-1.))
           dxdz =(size(1)*real(nxb)**(-1.)) * (size(3)*real(nzb)**(-1.)) 

           if (update_fluxes(lb)) then

              delx = nxb/size(1)
              dely = nyb/size(2)
              delz = nzb/size(3)
              do k = kli-2, kui+1
                 do j = jli-2, jui+1
                    do i = ili-2, iui+1
                       mgfluxx(i,j,k,lb) = delx * & 
                            &             ( 9./8.*(work(i+1,j,k,lb,1) - work(i,j,k,lb,1)) - &
                            &              1./24.*(work(i+2,j,k,lb,1) - work(i-1,j,k,lb,1)))
         
                       mgfluxy(i,j,k,lb) = dely * &
                            &             ( 9./8.*(work(i,j+1,k,lb,1) - work(i,j,k,lb,1)) - &
                            &              1./24.*(work(i,j+2,k,lb,1) - work(i,j-1,k,lb,1)))
 
                       mgfluxz(i,j,k,lb) = delz * &
                            &             ( 9./8.*(work(i,j,k+1,lb,1) - work(i,j,k,lb,1)) - &
                            &              1./24.*(work(i,j,k+2,lb,1) - work(i,j,k-1,lb,1)))
                    enddo
                 enddo
              enddo
           endif
           flux_x_ptr(FLXMC_FLUX,1,:,:,lb) = mgfluxx(ili-1,jli:jui,kli:kui,lb)*dydz
           flux_x_ptr(FLXMC_FLUX,2,:,:,lb) = mgfluxx(iui,jli:jui,kli:kui,lb)*dydz
           flux_y_ptr(FLXMC_FLUX,:,1,:,lb) = mgfluxy(ili:iui,jli-1,kli:kui,lb)*dxdz
           flux_y_ptr(FLXMC_FLUX,:,2,:,lb) = mgfluxy(ili:iui,jui,kli:kui,lb)*dxdz
           flux_z_ptr(FLXMC_FLUX,:,:,1,lb) = mgfluxz(ili:iui,jli:jui,kli-1,lb)*dxdy
           flux_z_ptr(FLXMC_FLUX,:,:,2,lb) = mgfluxz(ili:iui,jli:jui,kui,lb)*dxdy
        enddo
        

        end select



     endif

     !-------------------------------------------------------------------------------

     !               Now use the PARAMESH flux conservation routine to enforce
     !               continuity of the first derivative of the solution by
     !               overwriting the coarse-grid boundary flux at coarse-fine
     !               interfaces with the fine-grid boundary flux.

     if (flux_cons /= 0) then

        !               Call the PARAMESH flux conservation routine.
        call amr_flux_conserve(MyPE2, 0, 0)

        if (ndim .eq. 1) then
        do i = 1, lnblocks2
           mgfluxx(ili-1,jli:jui,kli:kui,i) = flux_x_ptr(FLXMC_FLUX,1,:,:,i)
           mgfluxx(iui,jli:jui,kli:kui,i)   = flux_x_ptr(FLXMC_FLUX,2,:,:,i)
        enddo
        elseif (ndim .eq. 2) then
           do i = 1, lnblocks2

              ! Get BlockSize:
              call Grid_getBlkPhysicalSize(i,size)

              dx = real(nxb)/size(1)
              dy = real(nyb)/size(2)

              mgfluxx(ili-1,jli:jui,kli:kui,i) = flux_x_ptr(FLXMC_FLUX,1,:,:,i)*dy
              mgfluxx(iui,jli:jui,kli:kui,i)   = flux_x_ptr(FLXMC_FLUX,2,:,:,i)*dy
              mgfluxy(ili:iui,jli-1,kli:kui,i) = flux_y_ptr(FLXMC_FLUX,:,1,:,i)*dx
              mgfluxy(ili:iui,jui,kli:kui,i)   = flux_y_ptr(FLXMC_FLUX,:,2,:,i)*dx
                  
           enddo

        elseif (ndim .eq. 3) then
           do i = 1, lnblocks2

              ! Get BlockSize:
              call Grid_getBlkPhysicalSize(i,size)

              dxdy =(real(nxb)*size(1)**(-1.)) * (real(nyb)*size(2)**(-1.))
              dydz =(real(nyb)*size(2)**(-1.)) * (real(nzb)*size(3)**(-1.))
              dxdz =(real(nxb)*size(1)**(-1.)) * (real(nzb)*size(3)**(-1.))


              mgfluxx(ili-1,jli:jui,kli:kui,i) = flux_x_ptr(FLXMC_FLUX,1,:,:,i)*dydz
              mgfluxx(iui,jli:jui,kli:kui,i)   = flux_x_ptr(FLXMC_FLUX,2,:,:,i)*dydz
              mgfluxy(ili:iui,jli-1,kli:kui,i) = flux_y_ptr(FLXMC_FLUX,:,1,:,i)*dxdz
              mgfluxy(ili:iui,jui,kli:kui,i)   = flux_y_ptr(FLXMC_FLUX,:,2,:,i)*dxdz
              mgfluxz(ili:iui,jli:jui,kli-1,i) = flux_z_ptr(FLXMC_FLUX,:,:,1,i)*dxdy
              mgfluxz(ili:iui,jli:jui,kui,i)   = flux_z_ptr(FLXMC_FLUX,:,:,2,i)*dxdy
           enddo
        endif

     endif

     !-------------------------------------------------------------------------------

     !               Now difference the fluxes and subtract from the RHS to
     !               obtain the residual.

     if (ndim == 1) then

        do lb = 1, lnblocks2
           if (update_values(lb)) then

              ! Get BlockSize:
              call Grid_getBlkPhysicalSize(lb,size)

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unkt,CENTER)      

              delx = nxb/size(1)
              do i = ili, iui
                 unkt(ires,i,jli:jui,kli:kui) = & 
                      &          unkt(irhs,i,jli:jui,kli:kui) - & 
                      &          delx * (mgfluxx(i,jli:jui,kli:kui,lb) - & 
                      &                  mgfluxx(i-1,jli:jui,kli:kui,lb))
              enddo

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unkt,CENTER)  

           endif
        enddo

     elseif (ndim == 2) then

        do lb = 1, lnblocks2
           if (update_values(lb)) then

              ! Get BlockSize:
              call Grid_getBlkPhysicalSize(lb,size)

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unkt,CENTER)      

              delx = real(nxb)/size(1)
              dely = real(nyb)/size(2)
              do j = jli, jui
                 do i = ili, iui
                    unkt(ires,i,j,kli:kui) = & 
                         &            unkt(irhs,i,j,kli:kui) - & 
                         &            delx * (mgfluxx(i,j,kli:kui,lb) - & 
                         &                    mgfluxx(i-1,j,kli:kui,lb)) - & 
                         &            dely * (mgfluxy(i,j,kli:kui,lb) - & 
                         &                    mgfluxy(i,j-1,kli:kui,lb))

                 !print*,"RES VALUE",lb,i,j,unkt(ires,i,j,kli:kui)

                 enddo
              enddo

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unkt,CENTER) 

           endif
        enddo

     else

        select case(gr_mgDiffOpDiscretize)

        case(2)  ! 2nd order Central difference Solution
        do lb = 1, lnblocks2
           if (update_values(lb)) then


              ! Get BlockSize:
              call Grid_getBlkPhysicalSize(lb,size)

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unkt,CENTER)  

              delx = nxb/size(1)
              dely = nyb/size(2)
              delz = nzb/size(3)
              do k = kli, kui
                 do j = jli, jui
                    do i = ili, iui
                       unkt(ires,i,j,k) = unkt(irhs,i,j,k) - & 
                            &              delx * (mgfluxx(i,j,k,lb) - mgfluxx(i-1,j,k,lb)) - & 
                            &              dely * (mgfluxy(i,j,k,lb) - mgfluxy(i,j-1,k,lb)) - & 
                            &              delz * (mgfluxz(i,j,k,lb) - mgfluxz(i,j,k-1,lb))
                    enddo
                 enddo
              enddo

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unkt,CENTER) 

           endif
        enddo

        case(4)    ! 4th order Central difference Solution
        do lb = 1, lnblocks2
           if (update_values(lb)) then

              ! Get BlockSize:
              call Grid_getBlkPhysicalSize(lb,size)

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unkt,CENTER)  

              delx = nxb/size(1)
              dely = nyb/size(2)
              delz = nzb/size(3)
              do k = kli, kui
                 do j = jli, jui
                    do i = ili, iui

                       unkt(ires,i,j,k) = unkt(irhs,i,j,k) - & 
                            &              delx * &
                            &             ( 9./8.*(mgfluxx(i,j,k,lb) - mgfluxx(i-1,j,k,lb))   - &
                            &              1./24.*(mgfluxx(i+1,j,k,lb) - mgfluxx(i-2,j,k,lb)))- &
                            &              dely * &
                            &             ( 9./8.*(mgfluxy(i,j,k,lb) - mgfluxy(i,j-1,k,lb))   - &
                            &              1./24.*(mgfluxy(i,j+1,k,lb) - mgfluxy(i,j-2,k,lb)))- &
                            &              delz * &
                            &             ( 9./8.*(mgfluxz(i,j,k,lb) - mgfluxz(i,j,k-1,lb))   - &
                            &              1./24.*(mgfluxz(i,j,k+1,lb) - mgfluxz(i,j,k-2,lb)))


                    enddo
                 enddo
              enddo

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unkt,CENTER) 

           endif
        enddo
        

        end select




     endif

     !===============================================================================

  end if

!  if (leaf_only.ne.0) then
!     call timer_stop("leaf_only")
!  else
!     call timer_stop("not leaf_only")
!  end if
!  call timer_stop("residual")
     
  return
end subroutine poisson_mg_residual



