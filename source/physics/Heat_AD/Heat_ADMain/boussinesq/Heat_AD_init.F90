subroutine Heat_AD_init(blockCount,blockList)

   use Heat_AD_data

   use IncompNS_data, only: ins_Ra, ins_Pr, &
                            ins_isCoupled, ins_isHeated, ins_meshME

   use RuntimeParameters_interface, only: RuntimeParameters_get, &
                                          RuntimeParameters_mapStrToInt
   implicit none

#include "constants.h"
#include "Flash.h"

   integer, INTENT(INOUT) :: blockCount
   integer, INTENT(INOUT) :: blockList(MAXBLOCKS)

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
   if (ins_isHeated) then
           ht_source_term = 1.0
   else
           ht_source_term = 0.0
   endif

   if (ins_meshMe .eq. MASTER_PE) then
     !write(*,*) 'Head_AD Boundary Conditions'
     !write(*,*) '            Low(Type / Value)                            High(Type / Value)'
     !write(*,*) 'I   : (',ht_Txl_type,' / ',ht_Txl_value,') (',ht_Txr_type,' / ',ht_Txr_value,')'
     !write(*,*) 'J   : (',ht_Tyl_type,' / ',ht_Tyl_value,') (',ht_Tyr_type,' / ',ht_Tyr_value,')'
     !write(*,*) 'K   : (',ht_Tzl_type,' / ',ht_Tzl_value,') (',ht_Tzr_type,' / ',ht_Tzr_value,')'
     !write(*,*) Txl_bcString,' ',Txr_bcString
     !write(*,*) Tyl_bcString,' ',Tyr_bcString
     !write(*,*) Tzl_bcString,' ',Tzr_bcString
     write(*,*) '1/SQRT(Ra Pr) = ',ht_invsqrtRaPr
     write(*,*) 'Volume Source = ',ht_source_term
   endif

end subroutine Heat_AD_init
