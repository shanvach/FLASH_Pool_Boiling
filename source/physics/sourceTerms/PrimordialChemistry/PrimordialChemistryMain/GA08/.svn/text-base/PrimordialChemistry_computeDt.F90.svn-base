!!****if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryMain/GA08/PrimordialChemistry_computeDt
!!
!!
!! NAME
!!
!!  PrimordialChemistry_computeDt
!!
!!
!! SYNOPSIS
!!
!!  call PrimordialChemistry_computeDt(integer(IN) :: blockID,
!!			   integer(IN) :: myPE,
!!			   real,pointer :: solnData(:,:,:,:),
!!			   real,(INOUT) :: dt_check,
!!			   integer(INOUT) :: dt_minloc(:) )
!!
!! DESCRIPTION
!!
!! Computes the timestep limiter
!!
!!
!! ARGUMENTS
!!
!! blockID	local block ID
!! myPE		local processor number
!! solnData	the physical, solution data from grid
!! dt_check	variable to hold timestep constraint
!! dt_minloc(5) array to hold limiting zone info: zone indices
!!
!!****

subroutine PrimordialChemistry_computeDt (blockID, myPE, blkLimits, blkLimitsGC, solnData, dt_chem, dt_minloc)

  use PrimordialChemistry_data
  use Simulation_data
  use pchem_networkInterface
  use pchem_interface
  use Grid_interface

#include "Flash.h"
#include "constants.h"

	implicit none

	integer, intent(IN) :: blockID, myPE
	integer, intent(IN), dimension(2,MDIM) :: blkLimits,blkLimitsGC
	real, INTENT(INOUT) :: dt_chem
	integer, INTENT(INOUT) :: dt_minloc(5)
	real, pointer :: solnData(:,:,:,:)
        integer		:: i,j,k,o
	real	:: dydx(NSPECIES), yy(NSPECIES), xx(NSPECIES)
	real    :: tlarge, ttstep, ttmin, ttt
	real    :: dt_temp, dt_check 
	real	:: temploc(5)
	integer :: specieMap
   tlarge = HUGE(0.0)

    do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
       do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
          do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
!	   ttt = 0.0
!	   !!Get species
!	   do o=1,NSPECIES
!              call pchem_mapNetworkToSpecies(o,specieMap)
!              xx(o) = solnData(specieMap,i,j,k)
!           enddo
!           call pchem_azbar()
!           do o=1,NSPECIES
!              yy(o) = ymass(o)
!           enddo
!           call pchem_network(ttt,yy,dydx)  !!should store odes into dydx
	   ttt = 1.0e20
!           !! check to see which is moving the fastest !!
!           do o=1,NSPECIES
!              ttmin = abs(sim_pchem_time*(yy(o) + 0.5*yy(iHP))/dydx(o))
           !   ttmin = 1.0e20 
!	     if(ttmin .lt. ttt) then
!                ttt = ttmin
!              endif
!           enddo
!
!	   call Grid_getBlkPtr(blockID,solnData)
!	   solnData(CHDT_VAR,i,j,k) = ttt
!	   call Grid_releaseBlkPtr(blockID,solnData)
!         
!           ttmin = 1.0e20
!
!           !!Check to make sure ttt is not less than zero
!           if(ttt .lt. 0.0) then
!		ttt = HUGE(0.0)
!	   endif
	   dt_check = ttt
!	   dt_temp = tlarge
!	   if(dt_check < dt_temp) then
		dt_temp = dt_check
		temploc(1) = i
		temploc(2) = j
		temploc(3) = k
		temploc(4) = blockID
		temploc(5) = MyPE
!	   endif
	   

	enddo
     enddo
   enddo


   if(dt_temp .lt. dt_chem) then
      dt_chem = dt_temp
      dt_minloc = temploc
   endif



	return
end subroutine PrimordialChemistry_computeDt
