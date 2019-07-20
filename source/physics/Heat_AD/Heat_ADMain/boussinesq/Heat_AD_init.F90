subroutine Heat_AD_init(blockCount,blockList)

   use Heat_AD_data

   use IncompNS_data, only: ins_Ra, ins_Pr, ins_isCoupled

   use RuntimeParameters_interface, only: RuntimeParameters_get, &
                                          RuntimeParameters_mapStrToInt
   implicit none

#include "constants.h"
#include "Flash.h"

   character(len=MAX_STRING_LENGTH) :: Txl_bcString, Txr_bcString
   character(len=MAX_STRING_LENGTH) :: Tyl_bcString, Tyr_bcString
   character(len=MAX_STRING_LENGTH) :: Tzl_bcString, Tzr_bcString

   ! Get the boundary conditions stored as strings 
   ! in the flash.par file
   call RuntimeParameters_get('Txl_boundary_type', Txl_bcString)
   call RuntimeParameters_get('Txr_boundary_type', Txr_bcString)
   call RuntimeParameters_get('Tyl_boundary_type', Tyl_bcString)
   call RuntimeParameters_get('Tyr_boundary_type', Tyr_bcString)
   call RuntimeParameters_get('Tzl_boundary_type', Tzl_bcString)
   call RuntimeParameters_get('Tzr_boundary_type', Tzr_bcString)

   ! Map the string boundary conditions to integer constants defined
   ! in the constants.h
   call RuntimeParameters_mapStrToInt(Txl_bcString, ht_Txl_type)
   call RuntimeParameters_mapStrToInt(Txr_bcString, ht_Txr_type)
   call RuntimeParameters_mapStrToInt(Tyl_bcString, ht_Tyl_type)
   call RuntimeParameters_mapStrToInt(Tyr_bcString, ht_Tyr_type)
   call RuntimeParameters_mapStrToInt(Tzl_bcString, ht_Tzl_type)
   call RuntimeParameters_mapStrToInt(Tzr_bcString, ht_Tzr_type)

   ! Get the boundary condition values stored in the flash.par file
   call RuntimeParameters_get('Txl_boundary_value', ht_Txl_value)
   call RuntimeParameters_get('Txr_boundary_value', ht_Txr_value)
   call RuntimeParameters_get('Tyl_boundary_value', ht_Tyl_value)
   call RuntimeParameters_get('Tyr_boundary_value', ht_Tyr_value)
   call RuntimeParameters_get('Tzl_boundary_value', ht_Tzl_value)
   call RuntimeParameters_get('Tzr_boundary_value', ht_Tzr_value)

   ht_invsqrtRaPr = 1. / (ins_Ra * ins_Pr)**0.5

end subroutine Heat_AD_init
