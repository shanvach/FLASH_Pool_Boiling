!!****if* source/Simulation/SimulationMain/magnetoHD/BB/hy_uhd_biermannSource
!!
!! NAME
!!
!!  hy_uhd_biermannSource
!!
!! SYNOPSIS
!!
!!  hy_uhd_biermannSource( integer (IN) :: blockCount,
!!                         integer (IN) :: blockList(blockCount),
!!                         real    (IN) :: dt )
!!
!! DESCRIPTION
!!
!! Implement Biermann Battery Term as a source to the magnetic field.
!!
!! ARGUMENTS
!!
!!  blockCount -  the number of blocks in blockList
!!  blockList  -  array holding local IDs of blocks on which to advance
!!  dt         -  timestep
!!***

Subroutine hy_uhd_biermannSource ( blockCount, blockList, dt )
  use Grid_interface, ONLY: Grid_getBlkIndexLimits, &
                            Grid_getBlkPtr,         &
                            Grid_releaseBlkPtr,     &
                            Grid_fillGuardCells

  use Hydro_data, ONLY : hy_biermannCoef,  &
                         hy_useBiermann,   &
                         hy_useBiermann1T, &
                         hy_bier1TZ,       &
                         hy_bier1TA,       &
                         hy_avogadro,      &
                         hy_qele,          &
                         hy_speedOfLight

  use Eos_interface, ONLY : Eos_wrapped,    &
                            Eos_getAbarZbar 
  
  implicit none

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  ! Arguments:
  integer, intent(IN) :: blockCount
  integer, intent(IN) :: blockList(blockCount)
  real,    intent(IN) :: dt

  ! Local Variables:
  real :: dens
  real :: abar
  real :: zbar
  real :: gradDensX, gradDensY
  real :: gradPresX, gradPresY
  real :: gradPeleX, gradPeleY
  real :: gradNeleX, gradNeleY
  real :: del(MDIM)
  real :: CL
  real :: QELE
  real :: NA
  real :: source
  real :: esource
  real :: ye
  real :: magGradNele, magGradPele
  real :: angle

  integer :: blockID
  integer :: blkLimits(LOW:HIGH,MDIM)
  integer :: blkLimitsGC(LOW:HIGH,MDIM)
  integer :: i, j, k, n, i1, j1

  real, pointer :: U(:,:,:,:)

  ! Begin:

  QELE = hy_qele
  CL = hy_speedOfLight
  NA = hy_avogadro


  do n = 1, blockCount
     blockID = blockList(n)

     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
     call Grid_getBlkPtr(blockID, U, CENTER)
     call Grid_getDeltas(blockID, del)

     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              ! Store nele everywhere:
              call Eos_getAbarZbar(U(:,i,j,k), Ye=ye)
              U(NELE_VAR,i,j,k) = ye * hy_avogadro * U(DENS_VAR,i,j,k)

              ! Store electric field everywhere:
              U(EX_VAR,i,j,k) = (U(PELE_VAR,i+1,j,k) - U(PELE_VAR,i-1,j,k))/(2*del(DIR_X)) / & 
                   (hy_qele * U(NELE_VAR,i,j,k))

              U(EY_VAR,i,j,k) = (U(PELE_VAR,i,j+1,k) - U(PELE_VAR,i,j-1,k))/(2*del(DIR_Y)) / & 
                   (hy_qele * U(NELE_VAR,i,j,k))
           end do
        end do
     end do

     call Grid_releaseBlkPtr(blockID,U,CENTER)
     call Eos_wrapped(MODE_DENS_EI_GATHER, blkLimits, blockID)
  end do

  call Grid_fillGuardCells(CENTER, ALLDIR)

  do n = 1, blockCount
     blockID = blockList(n)

     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
     call Grid_getBlkPtr(blockID, U, CENTER)
     call Grid_getDeltas(blockID, del)

     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              gradNelex = (U(NELE_VAR,i+1,j,k)-U(NELE_VAR,i-1,j,k))/(2.0*del(DIR_X))
              gradPeley = (U(PELE_VAR,i,j+1,k)-U(PELE_VAR,i,j-1,k))/(2.0*del(DIR_Y))
              
              gradNeley = (U(NELE_VAR,i,j+1,k)-U(NELE_VAR,i,j-1,k))/(2.0*del(DIR_Y))
              gradPelex = (U(PELE_VAR,i+1,j,k)-U(PELE_VAR,i-1,j,k))/(2.0*del(DIR_X))

              magGradNele = sqrt(gradNeleX**2 + gradNeleY**2)
              magGradPele = sqrt(gradPeleX**2 + gradPeleY**2)

              if(abs(magGradNele) < 1.0 .or. abs(magGradPele) < 1.0 .or. U(RHCD_VAR,i,j,k) < 0.5) then
                 U(ANGL_VAR,i,j,k) = 0.0
              else
                 U(ANGL_VAR,i,j,k) = 180.0/PI * acos( &
                      (gradPeleX*gradNeleX + gradPeleY*gradNeleY)/(magGradNele*magGradPele) )
              end if

           end do
        end do
     end do

     call Grid_releaseBlkPtr(blockID,U,CENTER)
     call Eos_wrapped(MODE_DENS_EI_GATHER, blkLimits, blockID)
  end do

  call Grid_fillGuardCells(CENTER, ALLDIR)

  do n = 1, blockCount
     blockID = blockList(n)

     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
     call Grid_getBlkPtr(blockID, U, CENTER)
     call Grid_getDeltas(blockID, del)

     call Eos_wrapped(MODE_DENS_EI_GATHER, blkLimits, blockID)

     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

#ifdef FLASH_3T
              ! do j1 = j-1,j+1
              !    do i1 = i-1,i+1
              !       ! Nele at i,j,k
              !       call Eos_getAbarZbar(U(:,i1,j1,k), Ye=ye)
              !       nele(2+i1-i,2+j1-j) = ye * hy_avogadro * U(DENS_VAR,i1,j1,k)
              !    end do
              ! end do
              
              gradNelex = (U(NELE_VAR,i+1,j,k)-U(NELE_VAR,i-1,j,k))/(2.0*del(DIR_X))
              gradPeley = (U(PELE_VAR,i,j+1,k)-U(PELE_VAR,i,j-1,k))/(2.0*del(DIR_Y))
              
              gradNeley = (U(NELE_VAR,i,j+1,k)-U(NELE_VAR,i,j-1,k))/(2.0*del(DIR_Y))
              gradPelex = (U(PELE_VAR,i+1,j,k)-U(PELE_VAR,i-1,j,k))/(2.0*del(DIR_X))

              magGradNele = sqrt(gradNeleX**2 + gradNeleY**2)
              magGradPele = sqrt(gradPeleX**2 + gradPeleY**2)

              ! ! Add Battery effect to Bz
              ! U(MAGZ_VAR,i,j,k) = U(MAGZ_VAR,i,j,k) &
              !      + dt*hy_biermannCoef/(hy_qele*nele(2,2)**2)*(gradPelex*gradNeley - gradPeley*gradNelex)

              ! ! Alternate method...
              ! ey_ip1_j = 1/(hy_qele*nele(3,2)*2*del(DIR_Y)) * (U(PELE_VAR,i+1,j+1,k)-U(PELE_VAR,i+1,j-1,k))
              ! ey_im1_j = 1/(hy_qele*nele(1,2)*2*del(DIR_Y)) * (U(PELE_VAR,i-1,j+1,k)-U(PELE_VAR,i-1,j-1,k))

              ! ex_i_jp1 = 1/(hy_qele*nele(2,3)*2*del(DIR_X)) * (U(PELE_VAR,i+1,j+1,k)-U(PELE_VAR,i-1,j+1,k))
              ! ex_i_jm1 = 1/(hy_qele*nele(2,1)*2*del(DIR_X)) * (U(PELE_VAR,i+1,j-1,k)-U(PELE_VAR,i-1,j-1,k))

              ! U(MAGZ_VAR,i,j,k) = U(MAGZ_VAR,i,j,k) + dt*hy_biermannCoef * &
              !        (ey_ip1_j - ey_im1_j) / (2*del(DIR_X)) &
              !      - (ex_i_jp1 - ex_i_jm1) / (2*del(DIR_Y))


              ! ! Alternate method...
              ! ey_ip1_j = 1/(hy_qele*nele(3,2)*2*del(DIR_Y)) * (U(PELE_VAR,i+1,j+1,k)-U(PELE_VAR,i+1,j-1,k))
              ! ey_im1_j = 1/(hy_qele*nele(1,2)*2*del(DIR_Y)) * (U(PELE_VAR,i-1,j+1,k)-U(PELE_VAR,i-1,j-1,k))

              ! ex_i_jp1 = 1/(hy_qele*nele(2,3)*2*del(DIR_X)) * (U(PELE_VAR,i+1,j+1,k)-U(PELE_VAR,i-1,j+1,k))
              ! ex_i_jm1 = 1/(hy_qele*nele(2,1)*2*del(DIR_X)) * (U(PELE_VAR,i+1,j-1,k)-U(PELE_VAR,i-1,j-1,k))

              ! if(U(RHCD_VAR,i,j,k) < 0.5 .and. i > 5 .and. j > 5) then
              U(BRMX_VAR,i,j,k) = hy_biermannCoef * &
                   magGradPele * magGradNele / (hy_qele * U(NELE_VAR,i,j,k)**2)
              
              ! U(ANGL_VAR,i,j,k) = 180.0/PI * acos( &
              !      (gradPeleX*gradNeleX + gradPeleY*gradNeleY)/(magGradNele*magGradPele) )
              
              if(U(RHCD_VAR,i,j,k) < 0.5) then
                 ! U(BRSC_VAR,i,j,k) = hy_biermannCoef * (  &
                 !      (U(EY_VAR,i+1,j,k) - U(EY_VAR,i-1,j,k)) / (2*del(DIR_X)) - &
                 !      (U(EX_VAR,i,j+1,k) - U(EX_VAR,i,j-1,k)) / (2*del(DIR_Y)) )

                 U(BRSC_VAR,i,j,k) = hy_biermannCoef / (hy_qele*U(NELE_VAR,i,j,k)**2) * &
                     (gradPelex*gradNeley - gradPeley*gradNelex)
                 

                 ! angle = U(ANGL_VAR,i,j,k) * PI/180.0

                 ! U(BRSC_VAR,i,j,k) = hy_biermannCoef / (hy_qele*U(NELE_VAR,i,j,k)**2) * &
                 !     magGradPele * magGradNele * sin(angle)

              else
                 U(BRSC_VAR,i,j,k) = 0.0
              end if

              
!     U(BRSC_VAR,i,j,k) = hy_biermannCoef / (hy_qele*U(NELE_VAR,i,j,k)**2) * &
!     (gradPelex*gradNeley - gradPeley*gradNelex)
              
!!              U(MAGZ_VAR,i,j,k) = U(MAGZ_VAR,i,j,k) + dt * U(BRSC_VAR,i,j,k)
              !! update the BB term only if the beta is high (>100) 
              if ( sqrt(2.*U(PRES_VAR,i,j,k)/(U(MAGZ_VAR,i,j,k) + dt*U(BRSC_VAR,i,j,k))**2) > 100 .and. U(BDRY_VAR,i,j,k)< 1.0) then            
!!              if (U(BDRY_VAR,i,j,k)< 1.0) then            
                U(MAGZ_VAR,i,j,k) = U(MAGZ_VAR,i,j,k) + dt * U(BRSC_VAR,i,j,k)  
              !! update the total energy to reflect the change in magnetic energy
                esource =dt*U(BRSC_VAR,i,j,k)*(2.*U(MAGZ_VAR,i,j,k)-dt*U(BRSC_VAR,i,j,k))*0.5/U(DENS_VAR,i,j,k)
                U(ENER_VAR,i,j,k) = U(ENER_VAR,i,j,k) - esource 
              endif

              
#endif
              ! if (hy_useBiermann1T) then

              !    dens = U(DENS_VAR,i,j,k)
              !    zbar = hy_bier1TZ
              !    abar = hy_bier1TA

              !    gradPresx=minmod(U(PRES_VAR,i+1,j,k)-U(PRES_VAR,i,j,  k),&
              !                     U(PRES_VAR,i,j,  k)-U(PRES_VAR,i-1,j,k))/del(DIR_X)

              !    gradPresy=minmod(U(PRES_VAR,i,j+1,k)-U(PRES_VAR,i,j,  k),&
              !                     U(PRES_VAR,i,j,  k)-U(PRES_VAR,i,j-1,k))/del(DIR_Y)

              !    gradDensx=minmod(U(DENS_VAR,i+1,j,k)-U(DENS_VAR,i,j,  k),&
              !                     U(DENS_VAR,i,j,  k)-U(DENS_VAR,i-1,j,k))/del(DIR_X)

              !    gradDensy=minmod(U(DENS_VAR,i,j+1,k)-U(DENS_VAR,i,j,  k),&
              !                     U(DENS_VAR,i,j,  k)-U(DENS_VAR,i,j-1,k))/del(DIR_Y)

              !    ! source stores the total change in the magnetic
              !    ! field for this cell:
              !    source = hy_biermannCoef * dt * &
              !         CL * abar / (QELE * NA * (1+zbar) *dens**2) * &
              !         (gradPresx * gradDensy - gradPresy * gradDensx)

              !    ! This source represents the conversion of internal
              !    ! to magnetic energy. There must be a corresponding
              !    ! decrease in einternal energy to account for this...

              !    ! E source is the specifc energy in ergs/g:
              !    esource = 0.5 * source**2 / U(DENS_VAR,i,j,k)                 
                 
              !    ! Update magnetic field and internal energy:
              !    U(MAGZ_VAR,i,j,k) = U(MAGZ_VAR,i,j,k) + source
              !    U(EINT_VAR,i,j,k) = U(EINT_VAR,i,j,k) - esource
              !    U(ENER_VAR,i,j,k) = U(ENER_VAR,i,j,k) - esource
              ! endif

           end do
        end do
     end do

     call Grid_releaseBlkPtr(blockID,U,CENTER)
     call Eos_wrapped(MODE_DENS_EI_GATHER, blkLimits, blockID)
  end do

  return

End Subroutine hy_uhd_biermannSource
