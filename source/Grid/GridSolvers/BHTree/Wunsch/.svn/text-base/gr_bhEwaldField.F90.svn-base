!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhEwaldField
!!
!! NAME
!!
!!  gr_bhEwaldField
!!
!!
!! SYNOPSIS
!!
!!   call gr_bhEwaldField()
!!
!! DESCRIPTION
!!
!!   Computes the contribution of a leaf block to the gravitational potential at a specific point.
!!
!! ARGUMENTS
!!
!!   none
!!
!!***




subroutine gr_bhEwaldField()
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  use gr_bhData, ONLY : gr_bhUseEwaldDecomp, gr_bhEwaldIsoFac, &
    gr_bhEwaldFieldNx, gr_bhEwaldFieldNy, gr_bhEwaldFieldNz, &
    gr_bhEwaldSeriesN, gr_bhEwaldFName, gr_bhEwaldAlwaysGenerate, &
    gr_bhLx, gr_bhLy, gr_bhLz, gr_bhEwaldXMax, gr_bhEwaldYMax, &
    gr_bhEwaldZMax, gr_bhTreeLrefineMax, gr_bhTreeMyPE, gr_bhComm, &
    gr_bhTreeNumProcs, gr_bhTreeEwald, gr_bhTreeCellSize

  implicit none
#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"

  real, parameter :: pi = PI
  integer :: i, j, k, ni, nj, nk, hi, hj, hk, ierr
  integer :: es_nx, es_ny, es_nz, i1d, chunk
  real :: xmax, xmin, ymax, ymin, zmax, zmin, Ln, Lav, dum
  real :: k2, kx, xni, yni, zni, rni, x, y, z, spii, ewald_alpha
  real :: Lx_per, Ly_per, Lz_per
  real :: ewald_1cell, re, rc, p
  real, allocatable :: loc_ewald(:,:,:)
  real :: et1, et2
  character(len=32) :: fname
  logical :: file_exists
  character(len=MAX_STRING_LENGTH) :: strBuff, grav_boundary_type_x &
  & , grav_boundary_type_y, grav_boundary_type_z

! If compilation fails in one of the next lines, or there are other problems
! related to ERFC, see file README.erfc in the BHTree/Wunsch implementation
! directory!

#ifdef USER_ERFC
#define ERFC USER_ERFC
  real, external :: USER_ERFC
#else
#define ERFC erfc
  intrinsic erfc
#endif

! If compilation fails in one of the preceding lines, or there are other problems
! related to ERFC, see file README.erfc in the BHTree/Wunsch implementation
! directory!


  call RuntimeParameters_get("grav_boundary_type_x", grav_boundary_type_x)
  call RuntimeParameters_get("grav_boundary_type_y", grav_boundary_type_y)
  call RuntimeParameters_get("grav_boundary_type_z", grav_boundary_type_z)
  call RuntimeParameters_get("gr_bhUseEwaldDecomp",  gr_bhUseEwaldDecomp)
  call RuntimeParameters_get("gr_bhEwaldIsoFac",     gr_bhEwaldIsoFac)

  ! return if periodic boundaries are not used in any direction
  if (    (grav_boundary_type_x .ne. "periodic") &
  & .and. (grav_boundary_type_y .ne. "periodic") &
  & .and. (grav_boundary_type_z .ne. "periodic")) return

!!$  if (gr_bhUseEwaldDecomp) then
!!$     ! Check here whether ERFC function is usable, abort if not??
!!$  end if

  call RuntimeParameters_get("gr_bhEwaldFieldNx", gr_bhEwaldFieldNx)
  call RuntimeParameters_get("gr_bhEwaldFieldNy", gr_bhEwaldFieldNy)
  call RuntimeParameters_get("gr_bhEwaldFieldNz", gr_bhEwaldFieldNz)
  call RuntimeParameters_get("gr_bhEwaldSeriesN", gr_bhEwaldSeriesN)
  call RuntimeParameters_get("gr_bhEwaldFName", gr_bhEwaldFName)
  call RuntimeParameters_get("gr_bhEwaldAlwaysGenerate", gr_bhEwaldAlwaysGenerate)
  call RuntimeParameters_get("xmin", xmin)
  call RuntimeParameters_get("xmax", xmax)
  call RuntimeParameters_get("ymin", ymin)
  call RuntimeParameters_get("ymax", ymax)
  call RuntimeParameters_get("zmin", zmin)
  call RuntimeParameters_get("zmax", zmax)
  allocate(gr_bhTreeEwald(0:gr_bhEwaldFieldNx-1, 0:gr_bhEwaldFieldNy-1, 0:gr_bhEwaldFieldNz-1))
  spii = 1.0 / sqrt(PI)
  gr_bhLx = xmax - xmin
  gr_bhLy = ymax - ymin
  gr_bhLz = zmax - zmin

  ! determine normalization constant Ln and average Lav
  if (grav_boundary_type_x == "periodic") then
    Lx_per = gr_bhLx
    es_nx = gr_bhEwaldSeriesN
    gr_bhEwaldXMax = 0.5 * gr_bhLx
  else
    Lx_per = gr_bhEwaldIsoFac*gr_bhLx
    es_nx = gr_bhEwaldSeriesN
    gr_bhEwaldXMax = gr_bhLx
  endif
  if (grav_boundary_type_y == "periodic") then
    Ly_per = gr_bhLy
    es_ny = gr_bhEwaldSeriesN
    gr_bhEwaldYMax = 0.5 * gr_bhLy
  else
    Ly_per = gr_bhEwaldIsoFac*gr_bhLy
    es_ny = gr_bhEwaldSeriesN
    gr_bhEwaldYMax = gr_bhLy
  endif
  if (grav_boundary_type_z == "periodic") then
    Lz_per = gr_bhLz
    es_nz = gr_bhEwaldSeriesN
    gr_bhEwaldZMax = 0.5 * gr_bhLz
  else
    Lz_per = gr_bhEwaldIsoFac*gr_bhLz
    es_nz = gr_bhEwaldSeriesN
    gr_bhEwaldZMax = gr_bhLz
  endif
  Lav = (Lx_per+Ly_per+Lz_per)/3.
  ewald_alpha = 2.0/Lav
  Ln = Lx_per*Ly_per*Lz_per

  !print *, "ef: ", gr_bhEwaldFieldNx, gr_bhEwaldFieldNy, gr_bhEwaldFieldNz
  !print *, "es: ", es_nx, es_ny, es_nz

  ! check the existence of the ewald_field file
  if (gr_bhTreeMyPE == MASTER_PE) inquire(file=gr_bhEwaldFName, exist=file_exists)
  call MPI_Bcast(file_exists, 1, MPI_LOGICAL, MASTER_PE, gr_bhComm, ierr)

  if ((file_exists) .and. (.not. gr_bhEwaldAlwaysGenerate)) then
    ! ewald_field file is present and there is no user request to regenerate it
    ! => read the ewald field from it
    if (gr_bhTreeMyPE == MASTER_PE) then
      open(unit=53, file=gr_bhEwaldFName, status='old')
      do k = 0,gr_bhEwaldFieldNz-1
        do j = 0,gr_bhEwaldFieldNy-1
          do i = 0,gr_bhEwaldFieldNx-1
            read(53,*) x, y, z, gr_bhTreeEwald(i,j,k), dum
          enddo
          read(53,*)
        enddo
        read(53,*)
      enddo
      close(unit=53)
    endif
    call MPI_Bcast(gr_bhTreeEwald, gr_bhEwaldFieldNx*gr_bhEwaldFieldNy*gr_bhEwaldFieldNz &
    &  , FLASH_REAL, MASTER_PE, gr_bhComm, ierr)

  else
    ! generate the ewald field
    allocate(loc_ewald(0:gr_bhEwaldFieldNx-1, 0:gr_bhEwaldFieldNy-1, 0:gr_bhEwaldFieldNz-1))
    if (gr_bhTreeMyPE == MASTER_PE) then
      write (strBuff, '("Computing Ewald correction field, box size = ", i4, "x", i4, "x", i4)') &
           gr_bhEwaldFieldNx, gr_bhEwaldFieldNy, gr_bhEwaldFieldNz
      call Logfile_stamp( strBuff, "[BHTree]")
    endif

    do k = 0,gr_bhEwaldFieldNz-1
      do j = 0,gr_bhEwaldFieldNy-1
        do i = 0,gr_bhEwaldFieldNx-1

          ! on all CPUs set each element to zero at first
          gr_bhTreeEwald(i,j,k) = 0.0
          loc_ewald(i,j,k) = 0.0

          ! calculate the 1D index: 0..ewald_field_z*gr_bhEwaldFieldNy*ewald_field_x-1
          i1d = k*gr_bhEwaldFieldNy*gr_bhEwaldFieldNx + j*gr_bhEwaldFieldNx + i
          chunk = max(1, gr_bhEwaldFieldNz*gr_bhEwaldFieldNy*gr_bhEwaldFieldNx &
                / gr_bhTreeNumProcs + 1)
          if ((i1d >= gr_bhTreeMyPE*chunk) .and. (i1d < (gr_bhTreeMyPE+1)*chunk)) then
            ! coordinates of the point
            x = i * gr_bhEwaldXMax / (gr_bhEwaldFieldNx-1)
            y = j * gr_bhEwaldYMax / (gr_bhEwaldFieldNy-1)
            z = k * gr_bhEwaldZMax / (gr_bhEwaldFieldNz-1)

            ! first term
            et1 = 0
            do ni = -es_nx,es_nx
              do nj = -es_ny,es_ny
                do nk = -es_nz,es_nz
                  xni = x + ni*Lx_per
                  yni = y + nj*Ly_per
                  zni = z + nk*Lz_per
                  rni = sqrt(xni*xni + yni*yni + zni*zni)
                  
                  if (gr_bhUseEwaldDecomp) then
                    loc_ewald(i,j,k) = loc_ewald(i,j,k) + ERFC(ewald_alpha*rni) / (rni + 1.0d-99)
                    et1 = et1 + ERFC(ewald_alpha*rni) / (rni + 1.0d-99)
                    !print *, "EF1: ", i,j,k, erfc(ewald_alpha*rni) / rni, ewald_alpha*rni, ewald_alpha, rni
                  else
                    loc_ewald(i,j,k) = loc_ewald(i,j,k) + 1.0 / (rni + 1.0d-99)
                  endif
                enddo
              enddo
            enddo
            
            ! second term
            if (gr_bhUseEwaldDecomp) then
              et2 = 0
              do hi = -es_nx,es_nx
                do hj = -es_ny,es_ny
                  do hk = -es_nz,es_nz
                    if ((hi .eq. 0) .and. (hj .eq. 0) .and. (hk .eq. 0)) cycle
                    k2 = 4*pi*pi*(hi*hi/Lx_per/Lx_per + hj*hj/Ly_per/Ly_per + hk*hk/Lz_per/Lz_per)
                    kx = 2*pi*(hi*x/Lx_per + hj*y/Ly_per + hk*z/Lz_per)
                    loc_ewald(i,j,k) = loc_ewald(i,j,k) + (4*pi/(k2*Ln+1d-99)) &
                    & * exp(-k2/(4*ewald_alpha*ewald_alpha)) * cos(kx)
                    !print *, "EF2: ", i,j,k, k2, kx, Ln, (4*pi/k2/Ln) &
                    !& , exp(-k2/(4*ewald_alpha*ewald_alpha)), cos(kx)
                    et2 = et2 + (4*pi/(k2*Ln+1d-99)) &
                    & * exp(-k2/(4*ewald_alpha*ewald_alpha)) * cos(kx)
 
                  enddo
                enddo
              enddo
              !print *, "EF: ", gr_bhTreeMyPE, i,j,k, et1, et2, 1.0/sqrt(x*x+y*y+z*z)
            endif
          endif
        enddo
      enddo
    enddo

    ! distribute chunks of the ewald filed array to all CPUs
    call MPI_AllReduce(loc_ewald, gr_bhTreeEwald, &
    gr_bhEwaldFieldNx*gr_bhEwaldFieldNy*gr_bhEwaldFieldNz &
    &  , FLASH_REAL, MPI_SUM, gr_bhComm, ierr)

    ! extrapolate surrounding values to center where formula yields 0/0
    x = gr_bhTreeCellSize(gr_bhTreeLrefineMax, IAXIS)
    y = gr_bhTreeCellSize(gr_bhTreeLrefineMax, JAXIS)
    z = gr_bhTreeCellSize(gr_bhTreeLrefineMax, KAXIS)
    ewald_1cell = 0.0
    do ni = -es_nx,es_nx
      do nj = -es_ny,es_ny
        do nk = -es_nz,es_nz
          xni = x + ni*Lx_per
          yni = y + nj*Ly_per
          zni = z + nk*Lz_per
          rni = sqrt(xni*xni + yni*yni + zni*zni)
          
          if (gr_bhUseEwaldDecomp) then
            ewald_1cell = ewald_1cell + ERFC(ewald_alpha*rni) / (rni + 1d-99)
          else
            ewald_1cell = ewald_1cell + 1.0 / (rni + 1d-99)
          endif
        enddo
      enddo
    enddo
    ! second term
    if (gr_bhUseEwaldDecomp) then
      do hi = -es_nx,es_nx
        do hj = -es_ny,es_ny
          do hk = -es_nz,es_nz
            if ((hi .eq. 0) .and. (hj .eq. 0) .and. (hk .eq. 0)) cycle
            k2 = 4*pi*pi*(hi*hi/Lx_per/Lx_per + hj*hj/Ly_per/Ly_per + hk*hk/Lz_per/Lz_per)
            kx = 2*pi*(hi*x/Lx_per + hj*y/Ly_per + hk*z/Lz_per)
            ewald_1cell = ewald_1cell + (4*pi/(k2*Ln+1d-99)) &
            & * exp(-k2/(4*ewald_alpha*ewald_alpha)) * cos(kx)
          enddo
        enddo
      enddo
    endif

    ! linear extrapolation
    re = sqrt((gr_bhEwaldXMax/(gr_bhEwaldFieldNx-1))**2 &
    &  +      (gr_bhEwaldYMax/(gr_bhEwaldFieldNy-1))**2 &
    &  +      (gr_bhEwaldZMax/(gr_bhEwaldFieldNz-1))**2)
    rc = sqrt(x*x+y*y+z*z)
    p = (-rc)/(re-rc)
    gr_bhTreeEwald(0,0,0) = (1.0-p)*ewald_1cell + p*gr_bhTreeEwald(1,1,1)

    ! write the ewald field into the file
    if (gr_bhTreeMyPE == MASTER_PE) then
      write (strBuff, '("Ewald correction field calculated.")')
      call Logfile_stamp( strBuff, "[BHTree]")
      
      open(unit=53, file=gr_bhEwaldFName, status='new')
      do k = 0,gr_bhEwaldFieldNz-1
        z = (k * 0.5 * gr_bhLz) / gr_bhEwaldFieldNz
        do j = 0,gr_bhEwaldFieldNy-1
          y = (j * 0.5 * gr_bhLy) / gr_bhEwaldFieldNy
          do i = 0,gr_bhEwaldFieldNx-1
            x = (i * 0.5 * gr_bhLx) / gr_bhEwaldFieldNx
            write(53,'(5(e17.10, 2x))') x, y, z, gr_bhTreeEwald(i,j,k) &
            , 1.0/(sqrt(x*x+y*y+z*z)+1d-99)
          enddo
          write(53,*)
        enddo
        write(53,*)
      enddo
      close(unit=53)
    endif ! MASTER_PE
    deallocate(loc_ewald)
  endif ! generate the ewald field

end subroutine gr_bhEwaldField

