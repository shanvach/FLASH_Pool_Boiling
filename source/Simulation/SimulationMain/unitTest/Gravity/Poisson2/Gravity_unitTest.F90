!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson2/Gravity_unitTest
!! NAME
!!
!!  Gravity_unitTest
!! 
!! SYNOPSIS
!!
!!  call Gravity_unitTest(integer(IN) :: fileUnit,
!!                    logical(OUT) :: perfect
!!
!! DESCRIPTION
!!
!! This function is the unit test for the Gravity unit. It is invoked in
!! the setup unitTest/Gravity. The Config file for Gravity unit test setup
!! requests a few extra variables in the main grid data structure for
!! Grid scope temporary storage. The Simulation_initBlock of the Gravity
!! unit test initializes density in the right place for the DENS_VAR
!! variable (see Flash.h for DENS_VAR, TEMP_VAR etc definitions), and
!! temperature and pressure in the extra storage space CTMP_VAR
!! and CPRS_VAR. The physical quantities at this point are not in
!! thermal equilibrium. 
!!
!! The Gravity_unit test starts by copying the initialized
!! temperature into the TEMP_VAR location and calling the
!! Gravity_wrapped function with eosMode = MODE_DENS_TEMP, where
!! density and temperature are given and pressure and energy are
!! calculated. Now PRES_VAR and EINT_VAR contain values of pressure
!! and internal energy that are in thermal equilibrium, and the pressure values
!! are not necessarily what was stored in the extra storage space
!! during intialization. 
!! 
!! At this point in time three quantities; temperature,
!! pressure and energy are saved in the extra storage requested by
!! the unitTest/Gravity setup, say OTMP_VAR, OPRS_VAR and OENT_VAR. Now
!! the Gravity_unitTest function calls Gravity_wrapped with eosMode =
!! MODE_DENS_PRES, followed by eosMode= MODE_DENS_EI.  If the
!! newly calculated values of temperature, pressure and energy are
!! the same as those saved in OTMP_VAR, OPRS_VAR and OENT_VAR, then
!! we can conclude that the Gravity is working in MODE_DENS_PRES and
!! MODE_DENS_EI modes. However, we still can't say anything about the
!! MODE_DENS_TEMP mode. So we repeat the process by copying CPRS_VAR
!! into PRES_VAR and calling Gravity_wrapped with MODE_DENS_PRES. We
!! again save the calculated values in the extra storage and make two
!! more Gravity_wrapped calls with the remaining two modes. This time if
!! the new and old values of variables compare, we can conclude that
!! MODE_DENS_TEMP works too, and hence the unit test is successful.
!!
!! A final test calculates the optional derivates by setting the mask
!! argument to Gravity equal to true.  The values of these arguments are not
!! tested in any way.  The unitTest simply makes sure that they can be calculated
!! without NaNs or the like.
!!
!!  ARGUMENTS 
!!   
!!   
!! 
!!   fileUnit : unit number for file opened by the unitTest/Gravity setup
!!              in which to write results of the test
!!
!!   perfect : indicates test ran without error is true.
!!
!!  PARAMETERS
!!
!!  eintSwitch  a rarely used switch which ensures that internal energy calculations 
!!        maintain sufficient precision. Important only if energyTotal is dominated 
!!        by energyKinetic.
!!
!!***

subroutine Gravity_unitTest( fileUnit, perfect)
  use Gravity_interface, ONLY : Gravity_potentialListOfBlocks
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPtr,&
                             Grid_releaseBlkPtr, Grid_getBlkIndexLimits

  implicit none

# include "constants.h"
# include "Flash.h"

  integer, intent(in) ::  fileUnit
  logical, intent(out) :: perfect
  integer :: localBlkCount, blockID
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer, dimension(MAXBLOCKS) :: blkList
  integer :: blkCount
  real :: potError,factorMin,factorMax,apotAbsMax,gpotAbsMax
  real, parameter :: orig_tolerance = 1e-9 !unused
  real, parameter :: tolerance = 0.01 !accept 1% differences
  integer :: ib,ie,jb,je,kb,ke,i

  real, pointer, dimension(:,:,:,:):: solnData

  call Grid_getListOfBlocks(LEAF,blkList,blkCount)

  call Gravity_potentialListOfBlocks(blkcount,blkList)

  apotAbsMax = TINY(1.0)  !1.0E-99
  gpotAbsMax = TINY(1.0)   !1.0E-99
  factorMin =  1.0E10
  factorMax = -1.0E10
  potError  =  0.0
  do i=1,blkCount
     blockID=blkList(i)
     call Grid_getBlkPtr(blockId,solnData)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     !! In Simulation_initBlock,
     !! temperature is initialized in CTMP_VAR and pressure is
     !! initialized in CPRS_VAR. We don't change these variables
     !! we copy them into the usual variable name as needed.

     ib=blkLimits(LOW,IAXIS)
     ie=blkLimits(HIGH,IAXIS)

     jb=blkLimits(LOW,JAXIS)
     je=blkLimits(HIGH,JAXIS)

     kb=blkLimits(LOW,KAXIS)
     ke=blkLimits(HIGH,KAXIS)

     potError = max(potError,maxval(abs(solnData(GPOT_VAR,ib:ie,jb:je,kb:ke)-&
                                        solnData(APOT_VAR,ib:ie,jb:je,kb:ke))))
     solnData(FACT_VAR,ib:ie,jb:je,kb:ke)=&
                          solnData(APOT_VAR,ib:ie,jb:je,kb:ke)/&
                          solnData(GPOT_VAR,ib:ie,jb:je,kb:ke)
     apotAbsMax = max(apotAbsMax,maxval(abs(solnData(APOT_VAR,ib:ie,jb:je,kb:ke))))
     gpotAbsMax = max(gpotAbsMax,maxval(abs(solnData(GPOT_VAR,ib:ie,jb:je,kb:ke))))
!     factor = max(factor,maxval(abs(solnData(FACT_VAR,ib:ie,jb:je,kb:ke))))
     factorMax = max(factorMax,maxval(solnData(FACT_VAR,ib:ie,jb:je,kb:ke)))
     factorMin = min(factorMin,minval(solnData(FACT_VAR,ib:ie,jb:je,kb:ke)))
     call Grid_releaseBlkPtr(blockId,solnData)
  end do

  print*,' max(|apot|),max(|gpot|):',apotAbsMax,gpotAbsMax

  if(potError < tolerance*min(apotAbsMax,gpotAbsMax) .and. &
       factorMin .GE. 1-tolerance .and. &
       factorMax .LE. 1+tolerance) then
     perfect=.true.
     print*,potError,factorMin,factorMax
     print*,'SUCCESS all tests were fine'
  else
     perfect=.false.
     print*,'FAILURE some tests failed',potError,factorMin,factorMax
     
  end if
  return
end subroutine Gravity_unitTest




