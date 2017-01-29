!*******************************************************************************

!  Routine:     gr_bicgMultiAx()

!  Description: Computes the product y=A*x


!  Parameters:  irhs        Right-hand side, index of x.
!               ilhs        Left-hand side, index of y.


subroutine gr_bicgMultiAx (irhs, ilhs, gcellflg)

  !===============================================================================

  use bicg_common, only : ili, iui, jli, jui, kli, kui, &
                          ile, iue, jle, jue, kle, kue

  use Grid_data, ONLY : gr_meshMe,gr_meshNumProcs

  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getBlkPtr,       &
                                Grid_releaseBlkPtr,   &
                                Grid_getBlkPhysicalSize, &
                                Grid_getListOfBlocks, &
                                Grid_getBlkIndexLimits

  use gr_bicgInterface, ONLY : gr_bicgBndry


  use physicaldata, only : flux_x,flux_y,flux_z
  use workspace, ONLY: work
  use tree, only : nodetype,lrefine
  use paramesh_dimensions, only : nguard_work
  use paramesh_interfaces, only :  amr_flux_conserve

  use Driver_data,    ONLY : dr_nstep


  implicit none
#include "constants.h"
#include "Flash.h"  
#include "mpif.h"

  integer, intent(in) :: irhs, ilhs
  logical, intent(in) :: gcellflg


  integer :: i, j, k, lb, ierr
  real    :: delx, dely, delz
!  real    :: norm
!  logical :: some_proc_needs_updating, my_proc_needs_updating
  logical, dimension(MAXBLOCKS) :: update_fluxes
!  logical, dimension(MAXBLOCKS) :: update_values

  integer, parameter       :: MAXDIMS2 = 3
!  integer, parameter       :: MFACES2 = 2*MAXDIMS2
  real, dimension(MDIM) :: size
!  real, dimension(MFACES2)  :: neigh2, neigh_type



  real    :: mgfluxx(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC, & 
       &                   MAXBLOCKS)
  real    :: mgfluxy(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC, & 
       &                   MAXBLOCKS)
  real    :: mgfluxz(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC, & 
       &                   MAXBLOCKS)

  real, pointer, save, dimension(:,:,:,:,:) :: flux_x_ptr, flux_y_ptr, flux_z_ptr

  real, pointer, save, dimension(:,:,:,:) :: unkt

  !- kpd -
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData,facevarx,facevary,facevarz
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

  integer blockcount,ii,jj,blockID

  real bsize(MDIM),coord(MDIM)

  integer blockList(MAXBLOCKS)

  !=============================================================================


  if (first_call) then
     first_call = .false.
     myPE2 = gr_meshMe
     num_PEs = gr_meshNumProcs
     flux_x_ptr => flux_x
     flux_y_ptr => flux_y
     flux_z_ptr => flux_z
  end if
  call Grid_getLocalNumBlks(lnblocks2)


  flux_x_ptr(FLXBC_FLUX,:,:,:,:) = 0.
  flux_y_ptr(FLXBC_FLUX,:,:,:,:) = 0.
  flux_z_ptr(FLXBC_FLUX,:,:,:,:) = 0.

!  call Grid_getListOfBlocks(LEAF,blockList,blockCount)    

  update_fluxes(:) = .false.
!  if (blockCount .gt. 0) update_fluxes(blockList(1:blockCount)) = .true.
  do lb = 1, lnblocks2     
    update_fluxes(lb) = (nodetype(lb) .eq. 1)
  enddo

 


! ***
!  call Timers_start("gr_bicgBndry")
  call gr_bicgBndry (irhs,nguard_work,gcellflg)
!  call Timers_stop("gr_bicgBndry")


  ! Set Value of first point,block,proc to zero
  !if (gr_meshMe .eq. 0)  work(ili,jli,kli,blockList(1),1) = 0.

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

        flux_x_ptr(FLXBC_FLUX,1,:,:,lb) = mgfluxx(ili-1,jli:jui,kli:kui,lb)
        flux_x_ptr(FLXBC_FLUX,2,:,:,lb) = mgfluxx(iui,jli:jui,kli:kui,lb)



     enddo

  elseif (ndim == 2) then


     do lb = 1, lnblocks2

        ! Get BlockSize:
        call Grid_getBlkPhysicalSize(lb,size)
        !- kpd - 
        call Grid_getBlkPtr(lb,facexData,FACEX)
        call Grid_getBlkPtr(lb,faceyData,FACEY)
        call Grid_getBlkPtr(lb,facevarx ,FACEX)
        call Grid_getBlkPtr(lb,facevary ,FACEY)


        dx = size(1)/real(nxb)
        dy = size(2)/real(nyb)

!if (gr_meshMe .eq. 0) print*,"Checkpoint A"

        if (update_fluxes(lb)) then    !- KPD - Only For Leaf Blocks

           !print*,"Block List inside poisson_mg_residual.F90. lb=",lb,"blockCount=",lnblocks2
           !print*,"poisson_mg_residualMG.F90 Flux Update Blocks:",gr_meshMe,level,lb

           delx = real(nxb)/size(1)
           dely = real(nyb)/size(2)
           do j = jli-1, jui
              do i = ili-1, iui

                 !- kpd - Added for variable density. Calculate face densities...
!                    if (lrefine(lb) == level .OR. &
!                        ((lrefine(lb) .lt. level)  .and. (nodetype_save(lb) == 1)) ) then
                       MdensXR = ( facexData(RH1F_FACE_VAR,i+1,j  ,1) + &
                                   facexData(RH2F_FACE_VAR,i+1,j  ,1) )   ! Inverse density on right face.
                       MdensYR = ( faceyData(RH1F_FACE_VAR,i  ,j+1,1) + &
                                   faceyData(RH2F_FACE_VAR,i  ,j+1,1) )   ! Inverse density on top face.
!                    else
!                 MdensXR = facevarx(idenvar,i+1,j  ,1)
!                 MdensYR = facevary(idenvar,i  ,j+1,1)
!                    end if

                 if (MdensXR.lt.1. .OR. MdensYR.lt.1. ) then
                    if (dr_nstep .gt. 1) then
                       print*,"ERROR 0 poisson_mg_residualMG.F90: Density Equals Zero at a FACE",lb,i,j,MdensXR,MdensYR
                    end if
                 elseif (MdensXR.gt.1000.0 .OR. MdensYR.gt.1000.0 ) then
                    if (dr_nstep .gt. 1) then
                       print*,"ERROR poisson_mg_residualMG.F90: Inverse Density Greater than 1.0",lb,i,j,MdensXR,MdensYR
                    end if
                 end if

                 !- kpd - Constant density implementation
                 !---------------------------------------
                 !mgfluxx(i,j,kli:kui,lb) = & 
                 !     &            delx * (work(i+1,j,kli:kui,lb,1) - & 
                 !     &                    work(i,j,kli:kui,lb,1))
                 !mgfluxy(i,j,kli:kui,lb) = & 
                 !     &            dely * (work(i,j+1,kli:kui,lb,1) - & 
                 !     &                    work(i,j,kli:kui,lb,1))
                 
                 !---------------------------------------------------------------------------
                 !- kpd - Variable density implementation
                 !        This is odd b/c mgflux is defined on right face of cell (not left),
                 !           this oddity is made up for (re-adjusted) later on.
                 !---------------------------------------------------------------------------
!                 mgfluxx(i,j,kli:kui,lb) = & 
!                      delx * (work(i+1,j,kli:kui,lb,1) - & 
!                                        work(i  ,j,kli:kui,lb,1))
!                      !MdensXR * delx * (work(i+1,j,kli:kui,lb,1) - & 
!                      !                  work(i  ,j,kli:kui,lb,1))
!                 mgfluxy(i,j,kli:kui,lb) = & 
!                      dely * (work(i,j+1,kli:kui,lb,1) - & 
!                                        work(i,j  ,kli:kui,lb,1))
!                      !MdensYR * dely * (work(i,j+1,kli:kui,lb,1) - & 
!                      !                  work(i,j  ,kli:kui,lb,1))

                    mgfluxx(i,j,kli:kui,lb) = &
                         &            MdensXR * delx * (work(i+1,j,kli:kui,lb,1) - &
                         &                    work(i,j,kli:kui,lb,1))
                    mgfluxy(i,j,kli:kui,lb) = &
                         &            MdensYR * dely * (work(i,j+1,kli:kui,lb,1) - &
                         &                    work(i,j,kli:kui,lb,1))


              enddo
           enddo

           flux_x_ptr(FLXBC_FLUX,1,:,:,lb) = mgfluxx(ili-1,jli:jui,kli:kui,lb)*dy
           flux_x_ptr(FLXBC_FLUX,2,:,:,lb) = mgfluxx(iui,jli:jui,kli:kui,lb)*dy
           flux_y_ptr(FLXBC_FLUX,:,1,:,lb) = mgfluxy(ili:iui,jli-1,kli:kui,lb)*dx
           flux_y_ptr(FLXBC_FLUX,:,2,:,lb) = mgfluxy(ili:iui,jui,kli:kui,lb)*dx

        endif

!        flux_x_ptr(FLXBC_FLUX,1,:,:,lb) = mgfluxx(ili-1,jli:jui,kli:kui,lb)*dy
!        flux_x_ptr(FLXBC_FLUX,2,:,:,lb) = mgfluxx(iui,jli:jui,kli:kui,lb)*dy
!        flux_y_ptr(FLXBC_FLUX,:,1,:,lb) = mgfluxy(ili:iui,jli-1,kli:kui,lb)*dx
!        flux_y_ptr(FLXBC_FLUX,:,2,:,lb) = mgfluxy(ili:iui,jui,kli:kui,lb)*dx

        !- kpd - 
        call Grid_releaseBlkPtr(lb,facexData,FACEX)
        call Grid_releaseBlkPtr(lb,faceyData,FACEY)
        call Grid_releaseBlkPtr(lb,facevarx ,FACEX)
        call Grid_releaseBlkPtr(lb,facevary ,FACEY)

     enddo

  else

!     select case(gr_biDiffOpDiscretize)

!        case(2)  ! 2nd order Central difference Solution
        do lb = 1, lnblocks2

           ! Get BlockSize:
           call Grid_getBlkPhysicalSize(lb,size)

           !- kpd - 
           call Grid_getBlkPtr(lb,facexData,FACEX)
           call Grid_getBlkPtr(lb,faceyData,FACEY)
           call Grid_getBlkPtr(lb,facezData,FACEZ)
           call Grid_getBlkPtr(lb,facevarx ,FACEX)
           call Grid_getBlkPtr(lb,facevary ,FACEY)
           call Grid_getBlkPtr(lb,facevarz ,FACEZ)

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
!                       if (lrefine(lb) == level .OR. &
!                            ((lrefine(lb) .lt. level)  .and. (nodetype_save(lb) == 1)) ) then
!
                          MdensXR = ( facexData(RH1F_FACE_VAR,i+1,j  ,k) + &
                                      facexData(RH2F_FACE_VAR,i+1,j  ,k) )   ! Inverse density on right face.
                          MdensYR = ( faceyData(RH1F_FACE_VAR,i  ,j+1,k) + &
                                      faceyData(RH2F_FACE_VAR,i  ,j+1,k) )   ! Inverse density on top face.
                          MdensZR = ( facezData(RH1F_FACE_VAR,i  ,j,k+1) + &
                                      facezData(RH2F_FACE_VAR,i  ,j,k+1) )   ! Inverse density on front face.
!                       else
!                          MdensXR = facevarx(idenvar,i+1,j  ,k)
!                          MdensYR = facevary(idenvar,i  ,j+1,k)
!                          MdensZR = facevarz(idenvar,i  ,j,k+1)
!                       end if

                       if (MdensXR.lt.1. .OR. MdensYR.lt.1. .OR. MdensZR.lt.1.) then
                          if (dr_nstep .gt. 1) then
                          print*,"ERROR 1 poisson_mg_residualMG.F90:",lb,i,j,k,MdensXR,MdensYR,MdensZR
                          end if
                       elseif (MdensXR.gt.1000.0 .OR. MdensYR.gt.1000.0 .OR. MdensZR.gt.1000.0) then
                          if (dr_nstep .gt. 1) then
                          print*,"ERROR poisson_mg_residualMG.F90: Inverse Density Greater than 1.0",lb,i,j,MdensXR,MdensYR
                          end if
                       end if

                       !!---------------------------------------------------------------------------
                       !!- kpd - Constant Density Implementation
                       !!---------------------------------------------------------------------------
                       !mgfluxx(i,j,k,lb) = & 
                       !     &              delx * (work(i+1,j,k,lb,1) - work(i,j,k,lb,1))
                       !mgfluxy(i,j,k,lb) = & 
                       !     &              dely * (work(i,j+1,k,lb,1) - work(i,j,k,lb,1))
                       !mgfluxz(i,j,k,lb) = & 
                       !     &              delz * (work(i,j,k+1,lb,1) - work(i,j,k,lb,1))

                       !---------------------------------------------------------------------------
                       !- kpd - Variable density implementation
                       !        This is odd b/c mgflux is defined on right face of cell (not left),
                       !           this oddity is made up for (re-adjusted) later on.
                       !---------------------------------------------------------------------------
                       mgfluxx(i,j,k,lb) = &
                                   MdensXR * delx * (work(i+1,j,k,lb,1) - &
                                                     work(i  ,j,k,lb,1))
                       mgfluxy(i,j,k,lb) = &
                                   MdensYR * dely * (work(i,j+1,k,lb,1) - &
                                                     work(i,j  ,k,lb,1))
                       mgfluxz(i,j,k,lb) = &
                                   MdensZR * delz * (work(i,j,k+1,lb,1) - &
                                                     work(i,j  ,k,lb,1))

                    enddo
                 enddo
              enddo
           endif
           flux_x_ptr(FLXBC_FLUX,1,:,:,lb) = mgfluxx(ili-1,jli:jui,kli:kui,lb)*dydz
           flux_x_ptr(FLXBC_FLUX,2,:,:,lb) = mgfluxx(iui,jli:jui,kli:kui,lb)*dydz
           flux_y_ptr(FLXBC_FLUX,:,1,:,lb) = mgfluxy(ili:iui,jli-1,kli:kui,lb)*dxdz
           flux_y_ptr(FLXBC_FLUX,:,2,:,lb) = mgfluxy(ili:iui,jui,kli:kui,lb)*dxdz
           flux_z_ptr(FLXBC_FLUX,:,:,1,lb) = mgfluxz(ili:iui,jli:jui,kli-1,lb)*dxdy
           flux_z_ptr(FLXBC_FLUX,:,:,2,lb) = mgfluxz(ili:iui,jli:jui,kui,lb)*dxdy

           !- kpd - 
           call Grid_releaseBlkPtr(lb,facexData,FACEX)
           call Grid_releaseBlkPtr(lb,faceyData,FACEY)
           call Grid_releaseBlkPtr(lb,facezData,FACEZ)
           call Grid_releaseBlkPtr(lb,facevarx ,FACEX)
           call Grid_releaseBlkPtr(lb,facevary ,FACEY)
           call Grid_releaseBlkPtr(lb,facevarz ,FACEZ)

        enddo

!!$        case(4)    ! 4th order Central difference Solution
!!$        do lb = 1, lnblocks2
!!$
!!$           ! Get BlockSize:
!!$           call Grid_getBlkPhysicalSize(lb,size)
!!$
!!$           dxdy =(size(1)*real(nxb)**(-1.)) * (size(2)*real(nyb)**(-1.))
!!$           dydz =(size(2)*real(nyb)**(-1.)) * (size(3)*real(nzb)**(-1.))
!!$           dxdz =(size(1)*real(nxb)**(-1.)) * (size(3)*real(nzb)**(-1.)) 
!!$
!!$           if (update_fluxes(lb)) then
!!$
!!$              delx = nxb/size(1)
!!$              dely = nyb/size(2)
!!$              delz = nzb/size(3)
!!$              do k = kli-2, kui+1
!!$                 do j = jli-2, jui+1
!!$                    do i = ili-2, iui+1
!!$                       mgfluxx(i,j,k,lb) = delx * & 
!!$                            &             ( 9./8.*(work(i+1,j,k,lb,1) - work(i,j,k,lb,1)) - &
!!$                            &              1./24.*(work(i+2,j,k,lb,1) - work(i-1,j,k,lb,1)))
!!$         
!!$                       mgfluxy(i,j,k,lb) = dely * &
!!$                            &             ( 9./8.*(work(i,j+1,k,lb,1) - work(i,j,k,lb,1)) - &
!!$                            &              1./24.*(work(i,j+2,k,lb,1) - work(i,j-1,k,lb,1)))
!!$ 
!!$                       mgfluxz(i,j,k,lb) = delz * &
!!$                            &             ( 9./8.*(work(i,j,k+1,lb,1) - work(i,j,k,lb,1)) - &
!!$                            &              1./24.*(work(i,j,k+2,lb,1) - work(i,j,k-1,lb,1)))
!!$                    enddo
!!$                 enddo
!!$              enddo
!!$           endif
!!$           flux_x_ptr(FLXBC_FLUX,1,:,:,lb) = mgfluxx(ili-1,jli:jui,kli:kui,lb)*dydz
!!$           flux_x_ptr(FLXBC_FLUX,2,:,:,lb) = mgfluxx(iui,jli:jui,kli:kui,lb)*dydz
!!$           flux_y_ptr(FLXBC_FLUX,:,1,:,lb) = mgfluxy(ili:iui,jli-1,kli:kui,lb)*dxdz
!!$           flux_y_ptr(FLXBC_FLUX,:,2,:,lb) = mgfluxy(ili:iui,jui,kli:kui,lb)*dxdz
!!$           flux_z_ptr(FLXBC_FLUX,:,:,1,lb) = mgfluxz(ili:iui,jli:jui,kli-1,lb)*dxdy
!!$           flux_z_ptr(FLXBC_FLUX,:,:,2,lb) = mgfluxz(ili:iui,jli:jui,kui,lb)*dxdy
!!$        enddo
!!$        
!!$
!!$     end select



  endif

!if (gr_meshMe .eq. 0) print*,"Checkpoint B"

     !-------------------------------------------------------------------------------

     !               Now use the PARAMESH flux conservation routine to enforce
     !               continuity of the first derivative of the solution by
     !               overwriting the coarse-grid boundary flux at coarse-fine
     !               interfaces with the fine-grid boundary flux.
!  call Timers_start("amr_flux_conserve")
  call amr_flux_conserve(MyPE2, 0, 0)
!  call Timers_stop("amr_flux_conserve")

  if (ndim .eq. 1) then
     do i = 1, lnblocks2
        mgfluxx(ili-1,jli:jui,kli:kui,i) = flux_x_ptr(FLXBC_FLUX,1,:,:,i)
        mgfluxx(iui,jli:jui,kli:kui,i)   = flux_x_ptr(FLXBC_FLUX,2,:,:,i)
     enddo
  elseif (ndim .eq. 2) then
     do i = 1, lnblocks2

        ! Get BlockSize:
        call Grid_getBlkPhysicalSize(i,size)

        dx = real(nxb)/size(1)
        dy = real(nyb)/size(2)

        mgfluxx(ili-1,jli:jui,kli:kui,i) = flux_x_ptr(FLXBC_FLUX,1,:,:,i)*dy
        mgfluxx(iui,jli:jui,kli:kui,i)   = flux_x_ptr(FLXBC_FLUX,2,:,:,i)*dy
        mgfluxy(ili:iui,jli-1,kli:kui,i) = flux_y_ptr(FLXBC_FLUX,:,1,:,i)*dx
        mgfluxy(ili:iui,jui,kli:kui,i)   = flux_y_ptr(FLXBC_FLUX,:,2,:,i)*dx
        
     enddo

  elseif (ndim .eq. 3) then
     do i = 1, lnblocks2

        ! Get BlockSize:
        call Grid_getBlkPhysicalSize(i,size)

        dxdy =(real(nxb)*size(1)**(-1.)) * (real(nyb)*size(2)**(-1.))
        dydz =(real(nyb)*size(2)**(-1.)) * (real(nzb)*size(3)**(-1.))
        dxdz =(real(nxb)*size(1)**(-1.)) * (real(nzb)*size(3)**(-1.))


        mgfluxx(ili-1,jli:jui,kli:kui,i) = flux_x_ptr(FLXBC_FLUX,1,:,:,i)*dydz
        mgfluxx(iui,jli:jui,kli:kui,i)   = flux_x_ptr(FLXBC_FLUX,2,:,:,i)*dydz
        mgfluxy(ili:iui,jli-1,kli:kui,i) = flux_y_ptr(FLXBC_FLUX,:,1,:,i)*dxdz
        mgfluxy(ili:iui,jui,kli:kui,i)   = flux_y_ptr(FLXBC_FLUX,:,2,:,i)*dxdz
        mgfluxz(ili:iui,jli:jui,kli-1,i) = flux_z_ptr(FLXBC_FLUX,:,:,1,i)*dxdy
        mgfluxz(ili:iui,jli:jui,kui,i)   = flux_z_ptr(FLXBC_FLUX,:,:,2,i)*dxdy
     enddo
  endif

  !-------------------------------------------------------------------------------
  !               Now difference the fluxes and subtract from the RHS to
  !               obtain the residual.

  if (ndim == 1) then

     do lb = 1, lnblocks2
        if (update_fluxes(lb)) then

              ! Get BlockSize:
              call Grid_getBlkPhysicalSize(lb,size)

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unkt,CENTER)      

              delx = nxb/size(1)
              do i = ili, iui
                 unkt(ilhs,i,jli:jui,kli:kui) = & 
                      &          delx * (mgfluxx(i,jli:jui,kli:kui,lb) - & 
                      &                  mgfluxx(i-1,jli:jui,kli:kui,lb))
              enddo

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unkt,CENTER)  

         endif
      enddo

   elseif (ndim == 2) then

      do lb = 1, lnblocks2
         if (update_fluxes(lb)) then

            ! Get BlockSize:
            call Grid_getBlkPhysicalSize(lb,size)

            ! Point to blocks center vars:
            call Grid_getBlkPtr(lb,unkt,CENTER)      

            delx = real(nxb)/size(1)
            dely = real(nyb)/size(2)
            do j = jli, jui
               do i = ili, iui
                    !unkt(ilhs,i,j,kli:kui) = & 
                    !     &            delx * (mgfluxx(i,j,kli:kui,lb) - & 
                    !     &                    mgfluxx(i-1,j,kli:kui,lb)) - & 
                    !     &            dely * (mgfluxy(i,j,kli:kui,lb) - & 
                    !     &                    mgfluxy(i,j-1,kli:kui,lb))
                    unkt(ilhs,i,j,kli:kui) = & 
                         &            delx * (mgfluxx(i,j,kli:kui,lb) - & 
                         &                    mgfluxx(i-1,j,kli:kui,lb)) + &   !KPD - Sign Change !!!!!!
                         &            dely * (mgfluxy(i,j,kli:kui,lb) - & 
                         &                    mgfluxy(i,j-1,kli:kui,lb))

               enddo
            enddo
            
            ! Point to blocks center vars:
            call Grid_releaseBlkPtr(lb,unkt,CENTER) 

         endif
      enddo
      
  
      !if (gr_meshMe .eq. 0) then
      !  lb = blockList(1)
      !  ! Point to blocks center vars:
      !  call Grid_getBlkPtr(lb,unkt,CENTER)
      !  unkt(ilhs,ili,jli,kli) = 0.
      !  call Grid_releaseBlkPtr(lb,unkt,CENTER)
      !endif   

   else

!        select case(gr_biDiffOpDiscretize)


!        case(2)  ! 2nd order Central difference Solution
        do lb = 1, lnblocks2
           if (update_fluxes(lb)) then


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
                       unkt(ilhs,i,j,k) =  & 
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

!!$        case(4)    ! 4th order Central difference Solution
!!$        do lb = 1, lnblocks2
!!$           if (update_fluxes(lb)) then
!!$
!!$              ! Get BlockSize:
!!$              call Grid_getBlkPhysicalSize(lb,size)
!!$
!!$              ! Point to blocks center vars:
!!$              call Grid_getBlkPtr(lb,unkt,CENTER)  
!!$
!!$              delx = nxb/size(1)
!!$              dely = nyb/size(2)
!!$              delz = nzb/size(3)
!!$              do k = kli, kui
!!$                 do j = jli, jui
!!$                    do i = ili, iui
!!$
!!$                       unkt(ilhs,i,j,k) =  & 
!!$                            &              delx * &
!!$                            &             ( 9./8.*(mgfluxx(i,j,k,lb) - mgfluxx(i-1,j,k,lb))   - &
!!$                            &              1./24.*(mgfluxx(i+1,j,k,lb) - mgfluxx(i-2,j,k,lb)))- &
!!$                            &              dely * &
!!$                            &             ( 9./8.*(mgfluxy(i,j,k,lb) - mgfluxy(i,j-1,k,lb))   - &
!!$                            &              1./24.*(mgfluxy(i,j+1,k,lb) - mgfluxy(i,j-2,k,lb)))- &
!!$                            &              delz * &
!!$                            &             ( 9./8.*(mgfluxz(i,j,k,lb) - mgfluxz(i,j,k-1,lb))   - &
!!$                            &              1./24.*(mgfluxz(i,j,k+1,lb) - mgfluxz(i,j,k-2,lb)))
!!$
!!$
!!$                    enddo
!!$                 enddo
!!$              enddo
!!$
!!$              ! Point to blocks center vars:
!!$              call Grid_releaseBlkPtr(lb,unkt,CENTER) 
!!$
!!$           endif
!!$        enddo
!!$        
!!$
!!$        end select

     endif

     !===============================================================================

  return
end subroutine gr_bicgMultiAx

