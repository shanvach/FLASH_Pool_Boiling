!!****if* source/Simulation/SimulationMain/unitTest/Burn/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_initBlock(integer(IN) :: blockID)
!!
!!
!!
!!
!! DESCRIPTION
!!
!!  Initialize the solution data for a single block.  Here we setup
!!  a block with varying temperature, density and composition -- 
!!  composition linearly varying from 1.00 of sim_compa at the -x end
!!  to 1.00 of compb at the +x end; density going from rho_min at
!!  the -y end to rho_max at the +y end; and temperature varying
!!  from sim_tempMin at the -z end to sim_tempMax at the +z end.   
!!
!!
!! ARGUMENTS
!!
!!  blockID -       the integer id of the block to initialize
!!
!!
!! PARAMETERS
!!
!!     smallx       smallest allowable abundance
!!     xmin         boundary of domain
!!     xmax         boundary of domain
!!     ymin         boundary of domain
!!     ymax         boundary of domain
!!     zmin         boundary of domain
!!     zmax         boundary of domain
!!     sim_tempMin        lowest temperature in domain
!!     sim_tempMax        highest temperature in domain
!!     rho_min      lowest density in domain
!!     rho_max      highest density in domain
!!     sim_compa        material at left of domain
!!     compb        material at right of domain
!!
!!***

subroutine Simulation_initBlock(blockId)

  use Simulation_data, only:  sim_iMin, sim_jMin, sim_kMin, sim_iMax, sim_jMax, sim_kMax,&
       sim_compIndexA, sim_compIndexB, &
       sim_smallX, sim_tempMin, sim_tempMax, sim_rhoMin, sim_rhoMax

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
    Grid_getCellCoords, Grid_putPointData, Grid_releaseBlkPtr

  use Eos_interface, ONLY:  Eos_wrapped

#include "constants.h"
#include "Flash.h"

  implicit none
  integer,intent(IN) :: blockID
            
  character(100)      :: sim_compa, sim_compb
  real                :: dpt, dpd, det, ded, c_v, c_p, pel, ne, eta
  real                :: abar, zbar
  real                :: gamma, xx, yy, zz
  real                        :: temp, dens, pres, ener
  integer             :: i, j, k, n
  real                :: r
  real                :: mu
  integer             :: sizeX,sizeY,sizeZ,istat
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC, eosRange
  integer, dimension(3) :: pos
  real :: entropy, dst, dsd

  real,allocatable,dimension(:) :: xCenter,yCenter,zCenter
  real :: p, gam, vx, vy, vz, e, game
  real,dimension(SPECIES_BEGIN:SPECIES_END) :: xn
  logical :: gcell = .true.

  real, pointer :: solnData(:,:,:,:)

  ! Get the block locations
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
  
  allocate(xCenter(sizeX),stat=istat)
  allocate(yCenter(sizeY),stat=istat)
  allocate(zCenter(sizeZ),stat=istat)
  
  if(NDIM>2) then
     call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell,zCenter,sizeZ)
  else
     zCenter=0.0
  end if

  if(NDIM>1) then
     call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell,yCenter,sizeY)
  else
     yCenter=0.0
  end if

  call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell,xCenter,sizeX)

!!  Get a pointer to the solution data
  call Grid_getBlkPtr(blockID,solnData)

  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     zz = zCenter(k) 

     if (NDIM >2) then
        dens = (zz-sim_kmin)/(sim_kmax-sim_kmin)*(sim_rhoMax-sim_rhoMin)+sim_rhoMin
     else
        dens = (sim_rhoMin+sim_rhoMax)/2.
     endif

     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        yy = yCenter(j)
        if (NDIM == 1) then
           temp = (sim_tempMin+sim_tempMax)/2.
        else
           !! Note that this formula produces negative temperatures in (some of) the guard cells
           !! Hence, only run the unitTest within the physical domain
           temp = (yy-sim_jmin)/(sim_jmax-sim_jmin)*(sim_tempMax-sim_tempMin)+sim_tempMin
        endif
        
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           xx = xCenter(i) 
                   
           pos(IAXIS)=i
           pos(JAXIS)=j
           pos(KAXIS)=k
           eosRange(LOW, IAXIS) = i   ! Eos is being done on a point by point basis
           eosRange(HIGH,IAXIS) = i
           eosRange(LOW, JAXIS) = j
           eosRange(HIGH,JAXIS) = j
           eosRange(LOW, KAXIS) = k
           eosRange(HIGH,KAXIS) = k

           ! There are only two species, ranging from SPECIES_BEGIN to SPECIES_END
           xn = sim_smallx   ! initialize to slightly nonzero
           xn(sim_compIndexA) = (xx-sim_imin)/(sim_imax-sim_imin)*(1.-3.*sim_smallx)
           xn(sim_compIndexA) = xn(sim_compIndexA) + sim_smallx
           
           xn(sim_compIndexB) = (1.-sim_smallx) - xn(sim_compIndexA)
           !! Fill up abundances for ALL species, not just A&B
           do n=SPECIES_BEGIN,SPECIES_END
              call Grid_putPointData(blockID,CENTER,n,EXTERIOR,pos,&
                                    xn(n))
           enddo
           
           ! Now fundamental quantities
           vx = 0.
           vy = 0.
           vz = 0.
           
           !! Fill up the grid with these values
           call Grid_putPointData(blockID, CENTER, VELX_VAR, EXTERIOR, pos, vx)
           call Grid_putPointData(blockID, CENTER, VELY_VAR, EXTERIOR, pos, vy)
           call Grid_putPointData(blockID, CENTER, VELZ_VAR, EXTERIOR, pos, vz)
           call Grid_putPointData(blockID, CENTER, DENS_VAR, EXTERIOR, pos, dens)
           call Grid_putPointData(blockID, CENTER, TEMP_VAR, EXTERIOR, pos, temp)

           !! Now call Eos with MOD_DENS_TEMP, which fills up solnData with PRES and EINT, GAMC
           call Eos_wrapped(MODE_DENS_TEMP,eosRange,blockID)

           !! calculate game and store
           game = solnData(PRES_VAR,i,j,k) / (solnData(EINT_VAR,i,j,k)*solnData(DENS_VAR,i,j,k)) + 1.0
           call Grid_putPointData(blockID, CENTER, GAME_VAR, EXTERIOR, pos, game)
                      
        enddo
     end do
  end do

! Cleanup
  call Grid_releaseBlkPtr(blockID,solnData)
  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)
  
  return
end subroutine Simulation_initBlock


