subroutine gr_pfftDataToBuf(messageNode, direction)

#include "constants.h"
#include "Pfft.h"
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_pfftMessageNode, ONLY : message_node
  use gr_pfftFragmentNode, ONLY : fragment_node
  use gr_pfftReconfigData, ONLY : pfft_oneRefLev, pfft_logMode, pfft_logUnit, &
       pfft_srcGridVar, pfft_solnGridVar
  use gr_pfftReconfigFn, ONLY : gr_pfftSrlGridToBuf, gr_pfftSrlPencilToBuf

  implicit none
  type(message_node), pointer  :: messageNode
  integer, intent(IN) :: direction
  type(fragment_node), pointer :: fragmentNode

  if (.not.(associated(messageNode))) then
     call Driver_abortFlash("[gr_pfftDataToBuf]: "//&
          "ERROR - Message node is not associated")
  end if


  !Eventually this code should be replaced with: 
  !map_impureFn(fnGridToBuf, voidStruct, fragmentList)

  fragmentNode => messageNode % fragmentList % H
  do while(associated(fragmentNode))


     if (direction == TO_PFFT) then
       
        if (fragmentNode % metadata % srcFlashBlockRefLev &
             == pfft_oneRefLev) then
        
           call gr_pfftSrlGridToBuf(fragmentNode, messageNode % buf)

        else if (fragmentNode % metadata % srcFlashBlockRefLev &
             > pfft_oneRefLev) then

           !RESTRICT
           call Driver_abortFlash("[gr_pfftDataToBuf]: "//&
                "ERROR: Restrict not yet coded")

        else if (fragmentNode % metadata % srcFlashBlockRefLev &
             < pfft_oneRefLev) then

           !PROLONG
           call Driver_abortFlash("[gr_pfftDataToBuf]: "//&
                "ERROR: Prolong not yet coded")
           
        end if

     else if (direction == FROM_PFFT) then

        if (fragmentNode % metadata % srcFlashBlockRefLev &
             == pfft_oneRefLev) then
        
           call gr_pfftSrlPencilToBuf(fragmentNode, messageNode % buf)

        else if (fragmentNode % metadata % srcFlashBlockRefLev &
             > pfft_oneRefLev) then

           !RESTRICT
           call Driver_abortFlash("[gr_pfftDataToBuf]: "//&
                "ERROR: Restrict not yet coded")

        else if (fragmentNode % metadata % srcFlashBlockRefLev &
             < pfft_oneRefLev) then

           !PROLONG
           call Driver_abortFlash("[gr_pfftDataToBuf]: "//&
                "ERROR: Prolong not yet coded")
           
        end if
        
     else
    
        call Driver_abortFlash("[gr_pfftDataToBuf]: "//&
             "ERROR: Pfft direction invalid")

     end if


     fragmentNode => fragmentNode % next
  end do

end subroutine gr_pfftDataToBuf
