!!****if* source/Simulation/SimulationMain/magnetoHD/Riemann2D/Simulation_initBlock
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
!!  This version sets up a 2D version of shock-tube type MHD Riemann problem.
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!
!! 
!!
!!***


subroutine Simulation_initBlock(blockId)

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
  real :: xx, yy, zz
  real, allocatable, dimension(:) :: xCoord,yCoord,zCoord
  real :: enerZone, ekinZone, eintZone
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData


  ! dump some output to stdout listing the paramters
!!$   if (sim_meshMe == MASTER_PE) then
!!$1    format (1X, 1P, 4(A7, E13.7, :, 1X))
!!$2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
!!$  endif

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord(sizeX), stat=istat)
  allocate(yCoord(sizeY), stat=istat)
  allocate(zCoord(sizeZ), stat=istat)

  xCoord  = 0.0
  yCoord  = 0.0
  zCoord  = 0.0

  if (NDIM == 3) call Grid_getCellCoords(KAXIS, blockID, CENTER,sim_gCell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords(JAXIS, blockID, CENTER,sim_gCell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockID, CENTER,     sim_gCell, xCoord,  sizeX)
!------------------------------------------------------------------------------

  call Grid_getBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
  endif
#endif

  ! Loop over cells in the block.
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     zz = zCoord(k)
     
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        yy = yCoord(j)

        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           xx  = xCoord (i)

           ! Multiple species
           !solnData(SPECIES_BEGIN,i,j,k)=1.0e0-(NSPECIES-1)*sim_smallX
           do n=SPECIES_BEGIN,SPECIES_END
              solnData(n,i,j,k)=sim_smallX
           enddo

           !! Region 1: upper right
           if (xx .gt. 0.5 .and. yy .gt. 0.5) then

              solnData(DENS_VAR,i,j,k)  = sim_dens1
              solnData(PRES_VAR,i,j,k)  = sim_pres1
              solnData(VELX_VAR,i,j,k)  = sim_velx1
              solnData(VELY_VAR,i,j,k)  = sim_vely1
              solnData(VELZ_VAR,i,j,k)  = sim_velz1
              solnData(MAGX_VAR,i,j,k)  = sim_magx1
              solnData(MAGY_VAR,i,j,k)  = sim_magy1
              solnData(MAGZ_VAR,i,j,k)  = sim_magz1
#if NFACE_VARS > 0
              if (sim_killdivb) then
                 facexData(MAG_FACE_VAR,i,j,k)= sim_magx1
                 faceyData(MAG_FACE_VAR,i,j,k)= sim_magy1
              endif
#endif

           !! Region 2: upper left
           elseif (xx .le. 0.5 .and. yy .gt. 0.5) then

              solnData(DENS_VAR,i,j,k)  = sim_dens2
              solnData(PRES_VAR,i,j,k)  = sim_pres2
              solnData(VELX_VAR,i,j,k)  = sim_velx2
              solnData(VELY_VAR,i,j,k)  = sim_vely2
              solnData(VELZ_VAR,i,j,k)  = sim_velz2
              solnData(MAGX_VAR,i,j,k)  = sim_magx2
              solnData(MAGY_VAR,i,j,k)  = sim_magy2
              solnData(MAGZ_VAR,i,j,k)  = sim_magz2
#if NFACE_VARS > 0
              if (sim_killdivb) then
                 facexData(MAG_FACE_VAR,i,j,k)= sim_magx2
                 faceyData(MAG_FACE_VAR,i,j,k)= sim_magy2
              endif
#endif
           !! Region 3: lower left
           elseif (xx .le. 0.5 .and. yy .le. 0.5) then

              solnData(DENS_VAR,i,j,k)  = sim_dens3
              solnData(PRES_VAR,i,j,k)  = sim_pres3
              solnData(VELX_VAR,i,j,k)  = sim_velx3
              solnData(VELY_VAR,i,j,k)  = sim_vely3
              solnData(VELZ_VAR,i,j,k)  = sim_velz3
              solnData(MAGX_VAR,i,j,k)  = sim_magx3
              solnData(MAGY_VAR,i,j,k)  = sim_magy3
              solnData(MAGZ_VAR,i,j,k)  = sim_magz3
#if NFACE_VARS > 0
              if (sim_killdivb) then
                 facexData(MAG_FACE_VAR,i,j,k)= sim_magx3
                 faceyData(MAG_FACE_VAR,i,j,k)= sim_magy3
              endif
#endif
           !! Region 4: lower right
           elseif (xx .gt. 0.5 .and. yy .le. 0.5) then

              solnData(DENS_VAR,i,j,k)  = sim_dens4
              solnData(PRES_VAR,i,j,k)  = sim_pres4
              solnData(VELX_VAR,i,j,k)  = sim_velx4
              solnData(VELY_VAR,i,j,k)  = sim_vely4
              solnData(VELZ_VAR,i,j,k)  = sim_velz4
              solnData(MAGX_VAR,i,j,k)  = sim_magx4
              solnData(MAGY_VAR,i,j,k)  = sim_magy4
              solnData(MAGZ_VAR,i,j,k)  = sim_magz4
#if NFACE_VARS > 0
              if (sim_killdivb) then
                 facexData(MAG_FACE_VAR,i,j,k)= sim_magx4
                 faceyData(MAG_FACE_VAR,i,j,k)= sim_magy4
              endif
#endif                  
           endif

           solnData(DIVB_VAR,i,j,k)  = 0.
           solnData(MAGP_VAR,i,j,k)  = 0.5*(solnData(MAGX_VAR,i,j,k) **2 &
                                          + solnData(MAGY_VAR,i,j,k) **2 &
                                          + solnData(MAGZ_VAR,i,j,k) **2)

           ! Compute the gas energy and set the gamma-values needed for the EOS
           ekinZone = 0.5 * dot_product(solnData(VELX_VAR:VELZ_VAR,i,j,k),&
                                        solnData(VELX_VAR:VELZ_VAR,i,j,k))

           ! specific internal energy
           eintZone = solnData(PRES_VAR,i,j,k)/(sim_gamma-1.)/solnData(DENS_VAR,i,j,k)

           ! total specific gas energy
           enerZone = eintZone + ekinZone

           ! Take a limit value
           enerZone = max(enerZone, sim_smallP)

           solnData(ENER_VAR,i,j,k)=enerZone
           solnData(EINT_VAR,i,j,k)=eintZone
           solnData(GAMC_VAR,i,j,k)=sim_gamma
           solnData(GAME_VAR,i,j,k)=sim_gamma

        enddo
     enddo
  enddo


  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
  endif
#endif


  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
  
end subroutine Simulation_initBlock
