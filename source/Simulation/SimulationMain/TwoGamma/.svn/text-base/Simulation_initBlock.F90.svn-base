!!****if* source/Simulation/SimulationMain/TwoGamma/Simulation_initBlock
!!
!! NAME
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!  call Simulation_initBlock(  integer(in) :: blockID)
!!                         
!!
!! DESCRIPTION   
!!     Initialize fluid properties (density, pressure, velocity, etc.) 
!!          in a single block 
!!          for the Two gamma problem. There are two fluids in the 
!!          problem that propagate together.
!!
!! ARGUMENTS
!!      blockID    -      the current block number to be filled
!!
!!
!!***


subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY : sim_xmin, sim_xmax, sim_ymin, sim_ymax, &
                              sim_small,sim_p0,sim_rho1,sim_rho2,&
                              sim_gam1,sim_gam2,sim_int1,sim_int2,&
                              sim_cvelx,sim_gammac1, sim_gammac2,sim_xpert,&
                              sim_temp1,sim_temp2
  use Grid_interface, ONLY : &
    Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_putPointData
  use Eos_interface, ONLY : Eos_wrapped


  implicit none
#include "constants.h"
#include "Flash.h"

  integer, intent(IN) :: blockID

  integer :: i, j, k, n     !loop counters
  
  real, allocatable, DIMENSION(:) :: xCenter

  real, pointer, dimension(:,:,:,:):: solnData

  integer, dimension(LOW:HIGH,MDIM)::blkLimits,blkLimitsGC
  integer :: sizeX
  logical :: gcell=.true.
  real :: massfr1, massfr2,kinEner,game1,game2
  integer,dimension(MDIM):: startingPos

  ! No need to initialize unk variable storage to zero here, Flash takes care of that.
  !

  ! set coordinates
  !
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  allocate(xCenter(sizeX))

  call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell,xCenter,sizeX)
  massfr1=sim_small
  massfr2=1.0e0-sim_small
  kinEner=0.5*(sim_cvelx**2)
  game1=sim_p0/sim_int1*sim_rho1 + 1.0e0
  game2=sim_p0/sim_int2*sim_rho2 + 1.0e0
! now fill the master arrays
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     startingPos(KAXIS)=k
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        startingPos(JAXIS)=j
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           startingPos(IAXIS)=i

           call Grid_putPointData(blockID,CENTER,PRES_VAR,EXTERIOR,&
                                  startingPos,sim_p0)
           call Grid_putPointData(blockID,CENTER,VELX_VAR,EXTERIOR,&
                                  startingPos,sim_cvelx)
!!           call Grid_putPointData(blockID,CENTER,VELY_VAR,EXTERIOR,&
!!                                  startingPos,0.0)
!!           call Grid_putPointData(blockID,CENTER,VELZ_VAR,EXTERIOR,&
!!                                  startingPos,0.0)

!  if less than the perturbation, then material one
           if(xCenter(i).lt.sim_xpert) then
              call Grid_putPointData(blockID,CENTER,EINT_VAR,EXTERIOR,&
                                  startingPos,sim_int1)
              call Grid_putPointData(blockID,CENTER,ENER_VAR,EXTERIOR,&
                                  startingPos,sim_int1 + kinEner)
              call Grid_putPointData(blockID,CENTER,DENS_VAR,EXTERIOR,&
                   startingPos,sim_rho1)
              call Grid_putPointData(blockID,CENTER,TEMP_VAR,EXTERIOR,&
                   startingPos,sim_temp1)
              call Grid_putPointData(blockID,CENTER,GAMC_VAR,EXTERIOR,&
                   startingPos,sim_gammac1)
              call Grid_putPointData(blockID,CENTER,GAME_VAR,EXTERIOR,&
                   startingPos,game1)
              call Grid_putPointData(blockID,CENTER,SPECIES_BEGIN,EXTERIOR,&
                   startingPos,massfr2)
              call Grid_putPointData(blockID,CENTER,SPECIES_BEGIN+1,EXTERIOR,&
                   startingPos,massfr1)

!  else greater than the perturbation, then material two

           else

              call Grid_putPointData(blockID,CENTER,EINT_VAR,EXTERIOR,&
                                  startingPos,sim_int2)
              call Grid_putPointData(blockID,CENTER,ENER_VAR,EXTERIOR,&
                                  startingPos,sim_int2 + kinEner)
              call Grid_putPointData(blockID,CENTER,DENS_VAR,EXTERIOR,&
                   startingPos,sim_rho2)
              call Grid_putPointData(blockID,CENTER,TEMP_VAR,EXTERIOR,&
                   startingPos,sim_temp2)
              call Grid_putPointData(blockID,CENTER,GAMC_VAR,EXTERIOR,&
                   startingPos,sim_gammac2)
              call Grid_putPointData(blockID,CENTER,GAME_VAR,EXTERIOR,&
                   startingPos,game2)
              call Grid_putPointData(blockID,CENTER,SPECIES_BEGIN,EXTERIOR,&
                   startingPos,massfr1)
              call Grid_putPointData(blockID,CENTER,SPECIES_BEGIN+1,EXTERIOR,&
                   startingPos,massfr2)

           endif
!
        enddo
     enddo
  enddo
  !     
  
#ifdef EELE_VAR
  ! get the integer index information for the current block
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  call Eos_wrapped(MODE_DENS_EI_SCATTER,blkLimits,blockId)
#endif

  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           startingPos(IAXIS) = i
           startingPos(JAXIS) = j
           startingPos(KAXIS) = k
#ifdef ERAD_VAR
           call Grid_putPointData(blockId, CENTER, ERAD_VAR, EXTERIOR, startingPos, 0.0  )   
#endif
#ifdef E3_VAR
           call Grid_putPointData(blockId, CENTER, E3_VAR,   EXTERIOR, startingPos, 0.0  )   
#endif

#ifdef PRAD_VAR
           call Grid_putPointData(blockId, CENTER, PRAD_VAR, EXTERIOR, startingPos, 0.0  )   
#endif
#ifdef TRAD_VAR
           call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, startingPos, 0.0  )   
#endif
        enddo
     enddo
  enddo

  deallocate(xCenter)

  return
end subroutine Simulation_initBlock
