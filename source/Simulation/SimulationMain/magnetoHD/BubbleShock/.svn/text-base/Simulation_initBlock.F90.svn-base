!!****if* source/Simulation/SimulationMain/magnetoHD/BubbleShock/Simulation_initBlock
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
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!!  Dai & Woodward, JCP, 142:331--369, 1998
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!  
!!
!!
!!
!!***

subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY : sim_gCell, sim_gam1,sim_gam2, &
                              sim_smallX, &
                              sim_killdivb,           &
                              sim_densT,sim_presT,sim_velxT,sim_velyT,&
                              sim_velzT,sim_magxT,sim_magyT,sim_magzT,&
                              sim_densB,sim_presB,sim_velxB,sim_velyB,&
                              sim_velzB,sim_magxB,sim_magyB,sim_magzB,&
                              sim_lposn,sim_bubbleRadius, sim_bubbleXCtr, sim_bubbleYCtr,   &
                              sim_bubbleZCtr,sim_bubbleDensity

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  integer, intent(in) :: blockID
  !!$ ---------------------------------

  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ
  real,allocatable,dimension(:) :: xCoord,yCoordL,yCoordR,yCoord,zCoord
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: enerZone, ekinZone, eintZone, taper, radius, r0
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
  real :: One,Zero

  ! dump some output to stdout listing the paramters
!!$   if (sim_meshMe == MASTER_PE) then
!!$1    format (1X, 1P, 4(A7, E13.7, :, 1X))
!!$2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
!!$  endif

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord (sizeX),stat=istat)
  allocate(yCoordL(sizeY),stat=istat)
  allocate(yCoordR(sizeY),stat=istat)
  allocate(yCoord (sizeY),stat=istat)
  allocate(zCoord (sizeZ),stat=istat)

  xCoord = 0.0
  yCoordL= 0.0
  yCoordR= 0.0
  yCoord = 0.0
  zCoord = 0.0

  call Grid_getCellCoords(IAXIS, blockID, CENTER,     sim_gCell, xCoord,  sizeX)
  call Grid_getCellCoords(JAXIS, blockID, LEFT_EDGE,  sim_gCell, yCoordL, sizeY)
  call Grid_getCellCoords(JAXIS, blockID, RIGHT_EDGE, sim_gCell, yCoordR, sizeY)
  call Grid_getCellCoords(JAXIS, blockID, CENTER,     sim_gCell, yCoord,  sizeY)
  call Grid_getCellCoords(KAXIS, blockID, CENTER,     sim_gCell, zCoord,  sizeZ)
  !------------------------------------------------------------------------------

  call Grid_getBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_getBlkPtr(blockID,facezData,FACEZ)
  endif
#endif

  ! Bubble is the fluid 1 (FLD1) and the ambient flow is the fluid 2 (FLD2)

  One  = 1.-sim_smallX
  Zero = sim_smallX 

  r0 = sim_bubbleRadius**.98
  ! Loop over cells in the block.
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

           ! Multiple species
           !solnData(SPECIES_BEGIN,i,j,k)=1.0e0-(NSPECIES-1)*sim_smallX
           !do n=SPECIES_BEGIN,SPECIES_END
           !   solnData(n,i,j,k)=sim_smallX
           !enddo

           ! Cell-centered values
           if (yCoordL(j) > sim_lposn) then
              ! Left states
              solnData(FLD1_SPEC,i,j,k) = Zero
              solnData(FLD2_SPEC,i,j,k) = One
              solnData(DENS_VAR,i,j,k)  = sim_densT
              solnData(VELX_VAR,i,j,k)  = sim_velxT
              solnData(VELY_VAR,i,j,k)  = sim_velyT
              solnData(VELZ_VAR,i,j,k)  = sim_velzT
              solnData(MAGX_VAR,i,j,k)  = sim_magxT
              solnData(MAGY_VAR,i,j,k)  = sim_magyT
              solnData(MAGZ_VAR,i,j,k)  = sim_magzT
              solnData(PRES_VAR,i,j,k)  = sim_presT
              solnData(GAMC_VAR,i,j,k)  = sim_gam2
              solnData(GAME_VAR,i,j,k)  = sim_gam2
              solnData(EINT_VAR,i,j,k)  = sim_presT/(sim_gam2-1.)/sim_densT

              ekinZone = 0.5 * dot_product(solnData(VELX_VAR:VELZ_VAR,i,j,k),&
                                           solnData(VELX_VAR:VELZ_VAR,i,j,k))
              solnData(ENER_VAR,i,j,k) = ekinZone + solnData(EINT_VAR,i,j,k)

           ! Cell face-centered variables for StaggeredMesh scheme
#if NFACE_VARS > 0
              if (sim_killdivb) then
                 facexData(MAG_FACE_VAR,i,j,k)= sim_magxT
                 faceyData(MAG_FACE_VAR,i,j,k)= sim_magyT
                 if (NDIM == 3) facezData(MAG_FACE_VAR,i,j,k)= sim_magzT
              endif
#endif

           elseif ((yCoordL(j) < sim_lposn) .and. (yCoordR(j) > sim_lposn)) then
              solnData(FLD1_SPEC,i,j,k) = Zero
              solnData(FLD2_SPEC,i,j,k) = One
              solnData(DENS_VAR,i,j,k)  =.5*(sim_densT + sim_densB)
              solnData(VELX_VAR,i,j,k)  =.5*(sim_velxT + sim_velxB)
              solnData(VELY_VAR,i,j,k)  =.5*(sim_velyT + sim_velyB)
              solnData(VELZ_VAR,i,j,k)  =.5*(sim_velzT + sim_velzB)
              solnData(MAGX_VAR,i,j,k)  =.5*(sim_magxT + sim_magxB)
              solnData(MAGY_VAR,i,j,k)  =.5*(sim_magyT + sim_magyB)
              solnData(MAGZ_VAR,i,j,k)  =.5*(sim_magzT + sim_magzB)
              solnData(PRES_VAR,i,j,k)  =.5*(sim_presT + sim_presB)
              solnData(GAMC_VAR,i,j,k)  = sim_gam2
              solnData(GAME_VAR,i,j,k)  = sim_gam2
              solnData(EINT_VAR,i,j,k)  = (sim_presT+sim_presB)/&
                                          (sim_gam2-1.)/(sim_densT+sim_densB)

              ekinZone = 0.5 * dot_product(solnData(VELX_VAR:VELZ_VAR,i,j,k),&
                                           solnData(VELX_VAR:VELZ_VAR,i,j,k))
              solnData(ENER_VAR,i,j,k) = ekinZone + solnData(EINT_VAR,i,j,k)

           ! Cell face-centered variables for StaggeredMesh scheme
#if NFACE_VARS > 0
!!$              if (sim_killdivb) then
!!$                 facexData(MAG_FACE_VAR,i,j,k) = sim_magxT
!!$                 if (xCoord(i) < sim_lposn) then
!!$                    faceyData(MAG_FACE_VAR,i,j,k) = sim_magyT
!!$                 elseif (xCoord(i) == sim_lposn) then
!!$                    faceyData(MAG_FACE_VAR,i,j,k) = .5*(sim_magyT+sim_magyB)
!!$                 else
!!$                    faceyData(MAG_FACE_VAR,i,j,k) = sim_magyB
!!$                 endif
!!$                 if (NDIM == 3) facezData(MAG_FACE_VAR,i,j,k) = sim_magzT
!!$              endif

              if (sim_killdivb) then
                 facexData(MAG_FACE_VAR,i,j,k) = sim_magxT
                 if (yCoord(j) > sim_lposn) then
                    faceyData(MAG_FACE_VAR,i,j,k) = sim_magyT
                 elseif (yCoord(j) == sim_lposn) then
                    faceyData(MAG_FACE_VAR,i,j,k) = .5*(sim_magyT+sim_magyB)
                 else
                    faceyData(MAG_FACE_VAR,i,j,k) = sim_magyB
                 endif
                 if (NDIM == 3) facezData(MAG_FACE_VAR,i,j,k) = sim_magzT
              endif
#endif

           else
              ! Right states
              solnData(FLD1_SPEC,i,j,k) = Zero
              solnData(FLD2_SPEC,i,j,k) = One
              solnData(DENS_VAR,i,j,k)  = sim_densB
              solnData(VELX_VAR,i,j,k)  = sim_velxB
              solnData(VELY_VAR,i,j,k)  = sim_velyB
              solnData(VELZ_VAR,i,j,k)  = sim_velzB
              solnData(MAGX_VAR,i,j,k)  = sim_magxB
              solnData(MAGY_VAR,i,j,k)  = sim_magyB
              solnData(MAGZ_VAR,i,j,k)  = sim_magzB
              solnData(PRES_VAR,i,j,k)  = sim_presB
              solnData(GAMC_VAR,i,j,k)  = sim_gam2
              solnData(GAME_VAR,i,j,k)  = sim_gam2
              solnData(EINT_VAR,i,j,k)  = sim_presB/(sim_gam2-1.)/sim_densB

              ekinZone = 0.5 * dot_product(solnData(VELX_VAR:VELZ_VAR,i,j,k),&
                                           solnData(VELX_VAR:VELZ_VAR,i,j,k))
              solnData(ENER_VAR,i,j,k) = ekinZone + solnData(EINT_VAR,i,j,k)

              !print*,'right pres=',solnData(PRES_VAR,i,j,k)

           ! Cell face-centered variables for StaggeredMesh scheme
#if NFACE_VARS > 0
              if (sim_killdivb) then
                 facexData(MAG_FACE_VAR,i,j,k) = sim_magxB
                 faceyData(MAG_FACE_VAR,i,j,k) = sim_magyB
                 if (NDIM == 3) facezData(MAG_FACE_VAR,i,j,k) = sim_magzB
              endif
#endif
           endif

           ! Bubble (Fluid 1) density
           if (NDIM == 2) then
              radius = sqrt( (xCoord(i)-sim_bubbleXCtr)**2 &
                            +(yCoord(j)-sim_bubbleYCtr)**2)
           elseif (NDIM == 3) then
              radius = sqrt( (xCoord(i)-sim_bubbleXCtr)**2 &
                            +(yCoord(j)-sim_bubbleYCtr)**2 &
                            +(zCoord(k)-sim_bubbleZCtr)**2)
           endif

           taper=(sim_bubbleRadius-radius)/(sim_bubbleRadius-r0)

           if (radius < r0) then
              ! the interior of the bubble is initialized with Fluid 1 property
              solnData(FLD1_SPEC,i,j,k) = One
              solnData(FLD2_SPEC,i,j,k) = Zero
              solnData(DENS_VAR,i,j,k)  = sim_bubbleDensity

              solnData(VELX_VAR,i,j,k)  = sim_velxB
              solnData(VELY_VAR,i,j,k)  = sim_velyB
              solnData(VELZ_VAR,i,j,k)  = sim_velzB
              solnData(MAGX_VAR,i,j,k)  = sim_magxB
              solnData(MAGY_VAR,i,j,k)  = sim_magyB
              solnData(MAGZ_VAR,i,j,k)  = sim_magzB
              solnData(PRES_VAR,i,j,k)  = sim_presB

              solnData(GAMC_VAR,i,j,k)  = sim_gam1
              solnData(GAME_VAR,i,j,k)  = sim_gam1
              solnData(EINT_VAR,i,j,k)  = sim_presB/(sim_gam1-1.)/sim_densB

              ekinZone = 0.5 * dot_product(solnData(VELX_VAR:VELZ_VAR,i,j,k),&
                                           solnData(VELX_VAR:VELZ_VAR,i,j,k))
              solnData(ENER_VAR,i,j,k) = ekinZone + solnData(EINT_VAR,i,j,k)

           elseif ((radius > r0).and.(radius < sim_bubbleRadius)) then
              solnData(FLD1_SPEC,i,j,k) = Zero + (One  - Zero)*taper
              solnData(FLD2_SPEC,i,j,k) = One  + (Zero - One )*taper
              solnData(DENS_VAR,i,j,k)  = sim_densB + (sim_bubbleDensity-sim_densB)*taper

              solnData(VELX_VAR,i,j,k)  = sim_velxB
              solnData(VELY_VAR,i,j,k)  = sim_velyB
              solnData(VELZ_VAR,i,j,k)  = sim_velzB
              solnData(MAGX_VAR,i,j,k)  = sim_magxB
              solnData(MAGY_VAR,i,j,k)  = sim_magyB
              solnData(MAGZ_VAR,i,j,k)  = sim_magzB
              solnData(PRES_VAR,i,j,k)  = sim_presB

              solnData(GAMC_VAR,i,j,k)  = sim_gam2
              solnData(GAME_VAR,i,j,k)  = sim_gam2
              solnData(EINT_VAR,i,j,k)  = sim_presB/(sim_gam2-1.)/sim_densB

              ekinZone = 0.5 * dot_product(solnData(VELX_VAR:VELZ_VAR,i,j,k),&
                                           solnData(VELX_VAR:VELZ_VAR,i,j,k))

              solnData(ENER_VAR,i,j,k) = ekinZone + solnData(EINT_VAR,i,j,k)
           endif

           solnData(MAGP_VAR,i,j,k) = .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                     solnData(MAGX_VAR:MAGZ_VAR,i,j,k))


           ! Compute the gas energy and set the gamma-values needed for the EOS
           !ekinZone = 0.5 * dot_product(solnData(VELX_VAR:VELZ_VAR,i,j,k),&
           !                             solnData(VELX_VAR:VELZ_VAR,i,j,k))

           ! specific internal energy
           !eintZone = solnData(PRES_VAR,i,j,k)/(sim_gamma-1.)/solnData(DENS_VAR,i,j,k)

           ! total specific gas energy
           !enerZone = eintZone + ekinZone

           ! Take a limit value
           !enerZone = max(enerZone, sim_smallP)

           !solnData(ENER_VAR,i,j,k)=enerZone
           !solnData(EINT_VAR,i,j,k)=eintZone
           !solnData(GAMC_VAR,i,j,k)=sim_gamma
           !solnData(GAME_VAR,i,j,k)=sim_gamma
           !print*,'i,j,x,y=',i,j,xCoord(i),yCoord(j)
           !print*,'gamc=',solnData(GAMC_VAR,i,j,k)
           !print*,''
        enddo
     enddo
  enddo

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
  endif
#endif

  deallocate(xCoord)
  deallocate(yCoordL)
  deallocate(yCoordR)
  deallocate(yCoord)
  deallocate(zCoord)

end subroutine Simulation_initBlock
