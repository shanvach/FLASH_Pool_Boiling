!!****if* source/Simulation/SimulationMain/DegenEOS/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  
!!  Initializes fluid data for a specified block.
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!
!! 
!!
!!***


subroutine Simulation_initBlock(blockID)

  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords, &
                             Grid_getBlkPtr,     &
                             Grid_releaseBlkPtr

  implicit none

#include "constants.h"
#include "Flash.h"


  !!$ Arguments -------------------
  integer, intent(in) :: blockID
  !!$------------------------------

  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: xx, xxL, xxR, yy, zz
  real, allocatable, dimension(:) :: xCoord,yCoord,zCoord, xCoordL, xCoordR
  real :: enerZone, ekinZone, eintZone
  real, pointer, dimension(:,:,:,:) :: solnData


  ! dump some output to stdout listing the paramters
  if (sim_meshMe == MASTER_PE) then
1    format (1X, 1P, 4(A7, E13.7, :, 1X))
2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
  endif

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord(sizeX), stat=istat)
  allocate(xCoordL(sizeX),stat=istat)  
  allocate(xCoordR(sizeX),stat=istat)
  allocate(yCoord(sizeY), stat=istat)
  allocate(zCoord(sizeZ), stat=istat)

  xCoord  = 0.0
  xCoordL = 0.0
  xCoordR = 0.0
  yCoord  = 0.0
  zCoord  = 0.0

  if (NDIM == 3) call Grid_getCellCoords(KAXIS, blockID, CENTER,sim_gCell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords(JAXIS, blockID, CENTER,sim_gCell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockID, CENTER,     sim_gCell, xCoord,  sizeX)
  call Grid_getCellCoords(IAXIS, blockID, LEFT_EDGE,  sim_gCell, xCoordL, sizeX)
  call Grid_getCellCoords(IAXIS, blockID, RIGHT_EDGE, sim_gCell, xCoordR, sizeX)
!------------------------------------------------------------------------------

  call Grid_getBlkPtr(blockID,solnData,CENTER)


  ! Loop over cells in the block.
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     !Cell left and right edges, width, and center (z-direction).
     
     zz = zCoord(k)
     
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        !Cell left and right edges, width, and center (y-direction).
        yy = yCoord(j)

        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           !Cell left and right edges, width, and center (x-direction).

!!$           ! Multiple species
!!$           solnData(SPECIES_BEGIN,i,j,k)=1.0e0-(NSPECIES-1)*sim_smallX
!!$           do n=SPECIES_BEGIN+1,SPECIES_END
!!$              solnData(n,i,j,k)=sim_smallX
!!$           enddo

           xx  = xCoord (i)
           xxL = xCoordL(i)
           xxR = xCoordR(i)

           if (xxR .lt. sim_posn) then
              ! Initialize cells to the left of the initial shock.

              solnData(DENS_VAR,i,j,k)  = sim_densH
              solnData(PRES_VAR,i,j,k)  = sim_pres0
              solnData(VELX_VAR,i,j,k)  = 0.
              solnData(VELY_VAR,i,j,k)  = 0.
              solnData(VELZ_VAR,i,j,k)  = 0.

              solnData(EINT_VAR,i,j,k)  = sim_eintH !1.607334e+18 !sim_eintH !1.607334e+18
              solnData(ENER_VAR,i,j,k)  = solnData(EINT_VAR,i,j,k)

              solnData(GAMC_VAR,i,j,k)  = sim_gamcH
              solnData(GAME_VAR,i,j,k)  = 1.+solnData(PRES_VAR,i,j,k)/&
                                          (solnData(DENS_VAR,i,j,k)*solnData(ENER_VAR,i,j,k ))

              solnData(TEMP_VAR,i,j,k)  = sim_tempH
              
           elseif ((xxR .ge. sim_posn) .and. (xxL .lt. sim_posn)) then
              ! Initialize cells which straddle the shock.  Treat them as
              ! though 1/2 of the cell lay to the left and 1/2 lay to the right.

              solnData(DENS_VAR,i,j,k)  = sim_densH/(1+sim_Atwood)
              solnData(PRES_VAR,i,j,k)  = sim_pres0
              solnData(VELX_VAR,i,j,k)  = 0.5*sqrt(sim_gamcL*solnData(PRES_VAR,i,j,k)/solnData(DENS_VAR,i,j,k))
              solnData(VELY_VAR,i,j,k)  = 0.
              solnData(VELZ_VAR,i,j,k)  = 0.

              solnData(EINT_VAR,i,j,k)  = 0.5*(sim_eintH+sim_eintL)
              solnData(ENER_VAR,i,j,k)  = solnData(EINT_VAR,i,j,k)

              solnData(GAMC_VAR,i,j,k)  = 0.5*(sim_gamcH+sim_gamcL)
              solnData(GAME_VAR,i,j,k)  = 1.+solnData(PRES_VAR,i,j,k)/&
                                          (solnData(DENS_VAR,i,j,k)*solnData(ENER_VAR,i,j,k ))

              solnData(TEMP_VAR,i,j,k)  = 0.5*(sim_tempH+sim_tempL)
              

              
           else
              ! Initialize cells to the right of the initial shock.

              solnData(DENS_VAR,i,j,k)  = sim_densH*(1.-sim_Atwood)/(1.+sim_Atwood)
              solnData(PRES_VAR,i,j,k)  = sim_pres0
              solnData(VELX_VAR,i,j,k)  = sqrt(sim_gamcL*solnData(PRES_VAR,i,j,k)/solnData(DENS_VAR,i,j,k))
              solnData(VELY_VAR,i,j,k)  = 0.
              solnData(VELZ_VAR,i,j,k)  = 0.

              solnData(EINT_VAR,i,j,k)  = sim_eintL !1.849076e+18 !sim_eintL !1.849076e+18
              solnData(ENER_VAR,i,j,k)  = solnData(EINT_VAR,i,j,k)

              solnData(GAMC_VAR,i,j,k)  = sim_gamcL
              solnData(GAME_VAR,i,j,k)  = 1.+solnData(PRES_VAR,i,j,k)/&
                                          (solnData(DENS_VAR,i,j,k)*solnData(ENER_VAR,i,j,k ))

              solnData(TEMP_VAR,i,j,k)  = sim_tempL


           endif

        enddo
     enddo
  enddo


  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)



  deallocate(xCoord)
  deallocate(xCoordL)
  deallocate(xCoordR)
  deallocate(yCoord)
  deallocate(zCoord)
  
end subroutine Simulation_initBlock
