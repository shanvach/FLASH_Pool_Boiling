module ib_viscoElastic_interface

        !implicit none

!INCLUDE .h files here
!#include "constants.h"
!#include "Flash.h"

        interface
        subroutine ib_dynamic_grid_advection(sn,sd,stest,adf,adfx,adfy,&
                                           cpt,dx,dy,dz,ix1,ix2,jy1,jy2)
        implicit none
        real, dimension(:,:,:), intent(inout) :: sn,adf,adfx,adfy
        real, dimension(:,:,:), intent(in)    :: sd,stest
        integer, dimension(:),  intent(in)    :: cpt
        real, intent(in)    :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2
        end subroutine ib_dynamic_grid_advection
        end interface

        interface
        subroutine ib_levelset_constantprojection(s,so,u,v,dx,dy,ix1,ix2,jy1,jy2)
        implicit none
        real, dimension(:,:,:), intent(inout) :: s
        real, dimension(:,:,:), intent(in)    :: so,u,v
        real, intent(in)    :: dx,dy
        integer, intent(in) :: ix1,ix2,jy1,jy2
        end subroutine ib_levelset_constantprojection
        end interface

        interface
        subroutine ib_levelset_linearprojection(s,so,sn,u,v,dx,dy,ix1,ix2,jy1,jy2)
        implicit none
        real, dimension(:,:,:), intent(inout) :: s
        real, dimension(:,:,:), intent(in)    :: u,v
        real, dimension(:,:,:), intent(in)    :: so,sn
        real, intent(in)    :: dx,dy
        integer, intent(in) :: ix1,ix2,jy1,jy2
        end subroutine ib_levelset_linearprojection
        end interface

        interface
        subroutine ib_solid_stress(sd,sY,sX,A,AT,A_inv,&
                                   Taux,Tauy,Taum,Tau,&
                                   dx,dy,ix1,ix2,jy1,jy2)
        implicit none
        real, dimension(:,:,:), intent(inout) :: Tau,Taux,Tauy,Taum
        real, dimension(:,:,:), intent(in)    :: sd,sY,sX,A,AT,A_inv
        real, intent(in)    :: dx,dy
        integer, intent(in) :: ix1,ix2,jy1,jy2
        end subroutine ib_solid_stress
        end interface

end module ib_viscoElastic_interface
