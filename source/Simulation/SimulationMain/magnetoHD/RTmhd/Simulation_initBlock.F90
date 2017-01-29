!!****if* source/Simulation/SimulationMain/magnetoHD/RTmhd/Simulation_initBlock
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

  use Simulation_data, ONLY : sim_gCell, sim_gamma,     &
                              sim_smallX, sim_smallP,   &
                              sim_killdivb,&
                              sim_Bx0, sim_By0, sim_Bz0,&
                              sim_epsilon, sim_xMin,    &
                              sim_xMax, sim_yMin, sim_yMax,&
                              sim_zMin, sim_zMax,&
                              sim_gconst, sim_rhoHeavy, sim_rhoLight,&
                              sim_number, sim_meshMe

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  integer, intent(in) :: blockID
  !!$ ---------------------------------

  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ
  integer :: nn, mm
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: enerZone, ekinZone, eintZone
  real :: peps,pwid,pwid2,amp,ymid
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
  logical, save :: FirstCall = .TRUE.
  real, allocatable, dimension(:,:) :: rand_A

  ! Make sure enough fluids have been defined.
  if (NSPECIES /=2) then
     call Driver_abortFlash("Simulation_initBlock needs two fluids!")
  endif

  ! dump some output to stdout listing the paramters
   if (sim_meshMe == MASTER_PE) then
1    format (1X, 1P, 4(A7, E13.7, :, 1X))
2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
  endif

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord(sizeX),stat=istat)
  allocate(yCoord(sizeY),stat=istat)
  allocate(zCoord(sizeZ),stat=istat)

  xCoord = 0.0
  yCoord = 0.0
  zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords(KAXIS,blockID,CENTER,sim_gCell,zCoord,sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords(JAXIS,blockID,CENTER,sim_gCell,yCoord,sizeY)
  call Grid_getCellCoords(IAXIS,blockID,CENTER,sim_gCell,xCoord,sizeX)
  !------------------------------------------------------------------------------

  call Grid_getBlkPtr(blockID,solnData,CENTER)

  ! Half of the domain height
  ymid = 0.5*(sim_yMin+sim_yMax)

  ! Background pressure, or pressure at top of box
  if (sim_number == 1) then
     ! FLASH's default setup
     peps = 1.0
  else
     ! Jim Stone's setup
     peps = 2.5
  endif

  ! Things for the velocity perturbation---
  ! (1) width of the velocity perturbation
  pwid = 0.5

  ! (2) y0 in Brucespeak
  pwid2 = pwid/2.0

  ! (3) rt_taper_z in Rayspeak
  nn   = 6

  ! (4) mode number
  mm   = 1

  if (sim_number == 3) then
     allocate(rand_A(sizeX,sizeY))
     call random_seed()
     call random_number(rand_A) ! random number over the ragnes [0,1]
     rand_A = rand_A - 0.5
     rand_A = rand_A*0.01
  endif

  ! Write a message to stdout describing the problem setup
  if (FirstCall) then
     FirstCall = .false.
     if (sim_meshMe .eq. MASTER_PE) then
        ! Check that the grid has the correct shape
        if (sim_xMax.ne.(pwid/2.0).or.sim_zMax.ne.(pwid/2.0)) then 
           write(6,*)'WARNING: incorrect grid shape for div V = 0'
           write(6,*)'xmax, pwid, pwid/4 =',sim_xMax,pwid,pwid/4.0
        endif
     endif
  endif



#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_getBlkPtr(blockID,facezData,FACEZ)
  endif
#endif

  ! Loop over cells in the block.
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

           ! Multiple species
           if (yCoord(j) < ymid) then
              solnData(HEAV_SPEC,i,j,k)=sim_smallX
              solnData(LIGH_SPEC,i,j,k)=1.-sim_smallX
              solnData(DENS_VAR,i,j,k) =sim_rhoLight

              if (sim_number == 1) then
                 ! FLASH's default setup
                 solnData(PRES_VAR,i,j,k) =sim_rhoLight*sim_gconst*(yCoord(j)-ymid)&
                                          +sim_rhoHeavy*sim_gconst*(ymid-sim_yMax) + peps
              else
                 ! Jim Stone's setup
                 solnData(PRES_VAR,i,j,k) =peps - 0.1*solnData(DENS_VAR,i,j,k)*yCoord(j)
              endif
           else
              solnData(HEAV_SPEC,i,j,k)=1.-sim_smallX
              solnData(LIGH_SPEC,i,j,k)=sim_smallX
              solnData(DENS_VAR,i,j,k) =sim_rhoHeavy

              if (sim_number == 1) then
                 ! FLASH's defaults setup
                 solnData(PRES_VAR,i,j,k) =sim_rhoHeavy*sim_gconst*(yCoord(j)-sim_yMax) + peps
              else
                 ! Jim Stone's setup
                 solnData(PRES_VAR,i,j,k) =peps - 0.1*solnData(DENS_VAR,i,j,k)*yCoord(j)
              endif
           endif

           solnData(TEMP_VAR,i,j,k) = 1.

           ! Cell-centered values
           ! velocity perturbation
           if (sim_number == 1) then
              if (abs(yCoord(j)-ymid) < pwid2/mm) then
                 amp = -sim_epsilon*sqrt(sim_gamma*solnData(PRES_VAR,i,j,k)/solnData(DENS_VAR,i,j,k))

#if NDIM == 2
                 solnData(VELX_VAR,i,j,k)= &
                      0.25 * amp * nn * sim_xMax / pwid2 *  &
                      sin(2.0*mm*PI*xCoord(i)/sim_xMax) * &
                      cos(0.5*mm*PI*(yCoord(j)-ymid+pwid2)/pwid2) *  &
                      (sin(0.5*mm*PI*(yCoord(j)-ymid+pwid2)/pwid2))**(nn - 1)

                 solnData(VELY_VAR,i,j,k)= &
                      -amp * &
                      cos(2.0*mm*PI*xCoord(i)/sim_xMax) * &
                      (sin(0.5*mm*PI*(yCoord(j)-ymid+pwid2)/pwid2))**nn

                 solnData(VELZ_VAR,i,j,k)= 0.

#elif NDIM == 3
                 solnData(VELX_VAR,i,j,k)= &
                      - 0.125 * amp * nn * sim_xMax / pwid2 *  &
                      sin(2.0*mm*PI*xCoord(i)/sim_xMax) * &
                      cos(2.0*mm*PI*zCoord(k)/sim_zMax) * &
                      cos(0.5*mm*PI*(yCoord(j)-ymid+pwid2)/pwid2) *  &
                      (sin(0.5*mm*PI*(yCoord(j)-ymid+pwid2)/pwid2))**(nn - 1)

                 solnData(VELY_VAR,i,j,k)=  &
                      amp * &
                      cos(2.0*mm*PI*xCoord(i)/sim_xMax) * &
                      cos(2.0*mm*PI*zCoord(k)/sim_zMax) * &
                      (sin(0.5*mm*PI*(yCoord(j)-ymid+pwid2)/pwid2))**nn

                 solnData(VELZ_VAR,i,j,k)=  &
                      - 0.125 * amp * nn * sim_zMax / pwid2 *  &
                      cos(2.0*mm*PI*xCoord(i)/sim_xMax) * &
                      sin(2.0*mm*PI*zCoord(k)/sim_zMax) * &
                      cos(0.5*mm*PI*(yCoord(j)-ymid+pwid2)/pwid2) *  &
                      (sin(0.5*mm*PI*(yCoord(j)-ymid+pwid2)/pwid2))**(nn - 1)

#endif
              else
                 solnData(VELX_VAR,i,j,k)=0.
                 solnData(VELY_VAR,i,j,k)=0.
                 solnData(VELZ_VAR,i,j,k)=0.
              endif
           else 
              ! Jim Stone's setup
              solnData(VELX_VAR,i,j,k) = 0.

              if (sim_number == 2) then
                 ! single-mode perturbation
                 solnData(VELY_VAR,i,j,k) = 0.01*(1.+cos(4.*PI*xCoord(i)))*(1.+cos(3.*PI*yCoord(j)))*0.25
              elseif (sim_number == 3) then
                 ! multi-mode perturbation
                 solnData(VELY_VAR,i,j,k) = rand_A(i,j)*(1.+cos(8./3.*PI*yCoord(j)))*0.5
              endif

              solnData(VELZ_VAR,i,j,k) = 0.
           endif

           solnData(GAMC_VAR,i,j,k)=sim_gamma
           solnData(GAME_VAR,i,j,k)=sim_gamma

#if defined(MAGX_VAR) && defined(MAGY_VAR) && defined(MAGZ_VAR) && defined(DIVB_VAR)
           solnData(MAGX_VAR,i,j,k)= sim_Bx0
           solnData(MAGY_VAR,i,j,k)= sim_By0
           solnData(MAGZ_VAR,i,j,k)= sim_Bz0
           solnData(MAGP_VAR,i,j,k)=  .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                     solnData(MAGX_VAR:MAGZ_VAR,i,j,k))
           solnData(DIVB_VAR,i,j,k)= 0.
#endif
           ! Compute the gas energy and set the gamma-values needed for the EOS
           ekinZone = 0.5 * dot_product(solnData(VELX_VAR:VELZ_VAR,i,j,k),&
                                        solnData(VELX_VAR:VELZ_VAR,i,j,k))

           ! specific internal energy
           eintZone = solnData(PRES_VAR,i,j,k)/(solnData(GAME_VAR,i,j,k)-1.)/solnData(DENS_VAR,i,j,k)

           ! total specific gas energy
           enerZone = eintZone + ekinZone

           ! Take a limit value
           enerZone = max(enerZone, sim_smallP)

           solnData(ENER_VAR,i,j,k)=enerZone
           solnData(EINT_VAR,i,j,k)=eintZone
           solnData(TEMP_VAR,i,j,k)=1.

        enddo
     enddo
  enddo



#if NFACE_VARS > 0
  !! Initialize face-centered magnetic fields
  if (sim_killdivb) then
     facexData(MAG_FACE_VAR,:,:,:) = sim_Bx0
     faceyData(MAG_FACE_VAR,:,:,:) = sim_By0
     if (NDIM==3) facezData(MAG_FACE_VAR,:,:,:) = sim_Bz0
  endif
#endif


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
  deallocate(yCoord)
  deallocate(zCoord)
  if (sim_number == 3) deallocate(rand_A)

end subroutine Simulation_initBlock
