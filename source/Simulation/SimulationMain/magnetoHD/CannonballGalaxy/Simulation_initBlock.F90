subroutine Simulation_initBlock(blockId)

  use Simulation_data

  use Eos_interface, ONLY : Eos

  use ut_interpolationInterface, ONLY : ut_hunt

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getDeltas, &
      Grid_putPointData, Grid_getCellCoords, Grid_putBlkData, &
      Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkRefineLevel
 
  implicit none

#include "Eos.h"
#include "constants.h"
#include "Flash.h"
  
  integer, INTENT(in) :: blockId

  ! Temporary variables.
   
  integer :: i, j, k, n, ii, jj, kk
  integer :: ir1, ir2, imin, imax, jmin, jmax, kmin, kmax
  integer :: blockRefineLevel

  real :: temp_zone, dens_zone, pres_zone, velx_zone
  real :: vely_zone, metl_zone, ener_zone, galx_zone
  real :: gamc_zone, velz_zone

  real :: sum_dens, sum_pres, sum_metl, sum_galx
  real :: pres_sample, dens_sample, metl_sample, Bmag, galx_sample
  real :: pres_sample1, dens_sample1, pres_center, dx_min
  real :: pres_sample2, dens_sample2, Btheta, Bphi, phi, theta
 
  real, save :: five_thirds = 5./3.
  real, save :: four_thirds = 4./3.
  real, save :: two_thirds = 2./3.

  real :: dr1, dr2, rho_aC, B_aC
  
  real, dimension(EOS_NUM) :: eosData
  logical,dimension(EOS_VARS+1:EOS_NUM) :: eosMask

  integer :: vecLen

  real, external :: interpolate

  ! Block coordinate info etc.

  real, dimension(MDIM)             :: del
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer, dimension(LOW:HIGH,MDIM) :: eosRange
  integer :: pos(3), size(3), ipos(3)

  real, dimension(:), allocatable :: xl, yl, zl, xc, yc, zc, xr, yr, zr

  real :: xx, yy, zz, dxx, dyy, dzz, dvol, rr, xx2, yy2, zz2, dx, dy, dz
  real, dimension(2) :: xgal, ygal, zgal, dgal
 
  real, dimension(:,:,:), allocatable :: rhog, t, vx, vy, vz, &
       galx, p, game, gamc, e, ei, bp, bx, by, bz, divb

  logical, save :: gcell = .true.

  real, external :: sim_getBfield
  real, pointer, dimension(:,:,:,:) :: facexData, faceyData, facezData

!===============================================================================

  pos(:) = 1

  eosMask = .false.

  call Grid_getBlkRefineLevel(blockID,blockRefineLevel)
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getDeltas(blockID,del)
  
  size(1) = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  size(2) = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  size(3) = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1

  allocate(rhog(size(1), size(2), size(3)))
  allocate(t(size(1), size(2), size(3)))
  allocate(vx(size(1), size(2), size(3)))
  allocate(vy(size(1), size(2), size(3)))
  allocate(vz(size(1), size(2), size(3)))
  allocate(galx(size(1), size(2), size(3)))
  allocate(e(size(1), size(2), size(3)))
  allocate(ei(size(1), size(2), size(3)))
  allocate(p(size(1), size(2), size(3)))
  allocate(game(size(1), size(2), size(3)))
  allocate(gamc(size(1), size(2), size(3)))
  allocate(bx(size(1), size(2), size(3)))
  allocate(by(size(1), size(2), size(3)))
  allocate(bz(size(1), size(2), size(3)))
  allocate(bp(size(1), size(2), size(3)))
  allocate(divb(size(1), size(2), size(3)))

  allocate(xl(size(1)), yl(size(2)), zl(size(3)))
  allocate(xc(size(1)), yc(size(2)), zc(size(3)))
  allocate(xr(size(1)), yr(size(2)), zr(size(3)))

  call Grid_getCellCoords(KAXIS, blockId, LEFT_EDGE, gcell, zl, size(3))
  call Grid_getCellCoords(JAXIS, blockId, LEFT_EDGE, gcell, yl, size(2))
  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, gcell, xl, size(1))
  call Grid_getCellCoords(KAXIS, blockId, CENTER, gcell, zc, size(3))
  call Grid_getCellCoords(JAXIS, blockId, CENTER, gcell, yc, size(2))
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xc, size(1))
  call Grid_getCellCoords(KAXIS, blockId, RIGHT_EDGE, gcell, zr, size(3))
  call Grid_getCellCoords(JAXIS, blockId, RIGHT_EDGE, gcell, yr, size(2))
  call Grid_getCellCoords(IAXIS, blockId, RIGHT_EDGE, gcell, xr, size(1))

#if NFACE_VARS > 0
  if (sim_killdivb .and. .not. sim_forceHydroLimit) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
  endif
#endif

  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)

     dzz = del(3) * nsubinv

     do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)

        dyy = del(2) * nsubinv
        
        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)

           dxx = del(1) * nsubinv
           
           sum_dens = 0.
           sum_pres = 0.
           sum_galx = 0.
           
           dvol = dxx*dyy*dzz
           
           do kk = 1, sim_subZones
              
              zz = zl(k) + (real(kk) - 0.5)*dzz
              zgal = zz - sim_zCtr
              
              do jj = 1, sim_subZones
                 
                 yy = yl(j) + (real(jj) - 0.5)*dyy
                 ygal = yy - sim_yCtr
           
                 do ii = 1, sim_subZones
                    
                    xx = xl(i) + (real(ii) - 0.5)*dxx
                    xgal = xx - sim_xCtr
                    
                    dgal = sqrt(xgal**2 + ygal**2 + zgal**2)
                    
                    pres_sample1 = 0.
                    pres_sample2 = 0.
                    dens_sample1 = 0.
                    dens_sample2 = 0.
                    
                    galx_sample  = 0.
               
                    if (.not. sim_testSingleGalaxy) then
                       
                       if (sim_ptdirn == 1) then
                          rr = xx
                       else if (sim_ptdirn == 2) then
                          rr = yy
                       else if (sim_ptdirn == 3) then
                          rr = zz
                       endif

                       dens_sample1 = interpolate(dens1, numPoints1, r1, rr)
                       pres_sample1 = interpolate(pres1, numPoints1, r1, rr)

                    endif

                    if (.not. sim_testAtmosphere) then
                       
                       dens_sample2 = interpolate(dens2, numPoints2, r2, dgal)
                       pres_sample2 = interpolate(pres2, numPoints2, r2, dgal)
                                              
                    endif

                    dens_sample = dens_sample1 + dens_sample2
                    pres_sample = pres_sample1 + pres_sample2

                    sum_dens = sum_dens + dens_sample
                    sum_pres = sum_pres + pres_sample
                    sum_galx = sum_galx + dens_sample2

                 enddo
              enddo
           enddo

           dens_zone = sum_dens*nsubvolinv
           pres_zone = sum_pres*nsubvolinv
           galx_zone = sum_galx/sum_dens

           velx_zone = 0.0
           vely_zone = 0.0
           velz_zone = 0.0
        
           if (.not. sim_testAtmosphere) then
              if (sim_ptdirn == 1) then
                 velx_zone = galx_zone*sim_vInit/dens_zone
              else if (sim_ptdirn == 2) then
                 vely_zone = galx_zone*sim_vInit/dens_zone                  
              else if (sim_ptdirn == 3) then
                 velz_zone = galx_zone*sim_vInit/dens_zone                 
              endif
           endif

           gamc_zone = five_thirds

           vx(i,j,k)   = velx_zone
           vy(i,j,k)   = vely_zone
           vz(i,j,k)   = velz_zone
           rhog(i,j,k) = dens_zone
           p(i,j,k)    = pres_zone
           game(i,j,k) = gamc_zone
           gamc(i,j,k) = gamc_zone

           eosData(EOS_PRES) = pres_zone
           eosData(EOS_DENS) = dens_zone

           vecLen = 1

           call Eos(MODE_DENS_PRES, vecLen, eosData)
           
           temp_zone   = eosData(EOS_TEMP)
           ener_zone   = eosData(EOS_EINT)

           t(i,j,k)    = temp_zone
           ei(i,j,k)   = ener_zone
           e(i,j,k)    = ener_zone + 0.5*(velx_zone**2 + &
                vely_zone**2+velz_zone**2)
           galx(i,j,k) = galx_zone

        enddo
     enddo
  enddo

#if NFACE_VARS > 0

    if (blockRefineLevel == sim_lrefineMin) then
  
      if (sim_killdivb .and. .not. sim_forceHydroLimit) then
     
         do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
                     
            do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        
               do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
                                                        
                  facexData(MAG_FACE_VAR,i,j,k) = sim_getBfield(1, xl(i), yl(j), zl(k))
                  faceyData(MAG_FACE_VAR,i,j,k) = sim_getBfield(2, xl(i), yl(j), zl(k))
                  facezData(MAG_FACE_VAR,i,j,k) = sim_getBfield(3, xl(i), yl(j), zl(k))

               enddo
            enddo
         enddo
     
      endif

   endif

   if (sim_killdivb .and. .not. sim_forceHydroLimit) then

      do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
         
         do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
            
            do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
               
               bx(i,j,k)   = 0.5*(facexData(MAG_FACE_VAR,i+1,j,k) + &
                    facexData(MAG_FACE_VAR,i,j,k))
               by(i,j,k)   = 0.5*(faceyData(MAG_FACE_VAR,i,j+1,k) + &
                    faceyData(MAG_FACE_VAR,i,j,k))
               bz(i,j,k)   = 0.5*(facezData(MAG_FACE_VAR,i,j,k+1) + &
                    facezData(MAG_FACE_VAR,i,j,k))
               
               divb(i,j,k) = (facexData(MAG_FACE_VAR,i+1,j,k) - &
                    facexData(MAG_FACE_VAR,i,j,k))/del(1) &
                    + (faceyData(MAG_FACE_VAR,i,j+1,k) - &
                    faceyData(MAG_FACE_VAR,i,j,k))/del(2) &
                    + (facezData(MAG_FACE_VAR,i,j,k+1) - &
                    facezData(MAG_FACE_VAR,i,j,k))/del(3)
              
#ifdef MAGP_VAR
               bp(i,j,k)   = 0.5*(bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2)
#endif

            enddo
         enddo
      enddo

   endif

#endif

#if NFACE_VARS > 0
  if (sim_killdivb .and. .not. sim_forceHydroLimit) then
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
  endif
#endif

  call Grid_putBlkData(blockID, CENTER, DENS_VAR, EXTERIOR, pos,rhog, size)
  call Grid_putBlkData(blockID, CENTER, PRES_VAR, EXTERIOR, pos,p,    size)
  call Grid_putBlkData(blockID, CENTER, GAME_VAR, EXTERIOR, pos,game, size)
  call Grid_putBlkData(blockID, CENTER, GAMC_VAR, EXTERIOR, pos,gamc, size)
  call Grid_putBlkData(blockID, CENTER, TEMP_VAR,     EXTERIOR, pos,t,    size)
  call Grid_putBlkData(blockID, CENTER, VELX_VAR,     EXTERIOR, pos,vx,   size)
  call Grid_putBlkData(blockID, CENTER, VELY_VAR,     EXTERIOR, pos,vy,   size)
  call Grid_putBlkData(blockID, CENTER, VELZ_VAR,     EXTERIOR, pos,vz,   size)
  call Grid_putBlkData(blockID, CENTER, GALX_MSCALAR, EXTERIOR, pos,galx, size)
  call Grid_putBlkData(blockID, CENTER, ENER_VAR,     EXTERIOR, pos,e,    size)
  call Grid_putBlkData(blockID, CENTER, EINT_VAR,     EXTERIOR, pos,ei,   size)
  call Grid_putBlkData(blockID, CENTER, MAGP_VAR,     EXTERIOR, pos,bp,   size)
  call Grid_putBlkData(blockID, CENTER, MAGX_VAR, EXTERIOR, pos,bx,   size)
  call Grid_putBlkData(blockID, CENTER, MAGY_VAR, EXTERIOR, pos,by,   size)
  call Grid_putBlkData(blockID, CENTER, MAGZ_VAR, EXTERIOR, pos,bz,   size)
  call Grid_putBlkData(blockID, CENTER, DIVB_VAR, EXTERIOR, pos,divb,   size)

  deallocate(rhog)
  deallocate(t)
  deallocate(vx)
  deallocate(vy)
  deallocate(vz)
  deallocate(galx)
  deallocate(e)
  deallocate(ei)
  deallocate(p)
  deallocate(gamc)
  deallocate(game)
  deallocate(bp)
  deallocate(bx)
  deallocate(by)
  deallocate(bz)
  deallocate(divb)
      
  deallocate(xl)
  deallocate(yl)
  deallocate(zl)
  deallocate(xc)
  deallocate(yc)
  deallocate(zc)
  deallocate(xr)
  deallocate(yr)
  deallocate(zr)
  
  return

end subroutine Simulation_initBlock

function sim_getBfield(dir, xx, yy, zz)

  use Simulation_data
  use ut_interpolationInterface, ONLY: ut_hunt

  implicit none

  integer, intent(IN) :: dir
  real, intent(IN) :: xx
  real, intent(IN) :: yy
  real, intent(IN) :: zz

  real :: sim_getBfield, Bfield, xi, yj, zk

  integer :: ii, jj, kk
     
  call ut_hunt(sim_Bxcoord, sim_nBzones, xx, ii)
  call ut_hunt(sim_Bycoord, sim_nBzones, yy, jj)
  call ut_hunt(sim_Bzcoord, sim_nBzones, zz, kk)

  if (ii < 1) ii = 1
  if (jj < 1) jj = 1
  if (kk < 1) kk = 1
  if (ii > sim_nBzones) ii = sim_nBzones
  if (jj > sim_nBzones) jj = sim_nBzones
  if (kk > sim_nBzones) kk = sim_nBzones

  if (dir == 1) then
     Bfield = sim_Bx(ii,jj,kk)
  else if (dir == 2) then
     Bfield = sim_By(ii,jj,kk)
  else if (dir == 3) then
     Bfield = sim_Bz(ii,jj,kk)
  endif

  sim_getBfield = Bfield

  return

end function sim_getBfield
