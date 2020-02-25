module ib_viscoElastic_interface

        interface
        subroutine ib_dynamic_grid_advection(sd,stest,&
                                             ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
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
        subroutine ib_solid_stress(sd,sX,sY,Tau1,Tau2,Tau3,Tau4&
                                   ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
        implicit none
        real, dimension(:,:,:), intent(inout) :: Tau1, Tau2, Tau3, Tau4
        real, dimension(:,:,:), intent(in)    :: sd,sX,sY
        real, dimension(:,:,:), intent(in)    :: A1,A2,A3,A4
        real, dimension(:,:,:), intent(in)    :: AT1,AT2,AT3,AT4
        real, dimension(:,:,:), intent(in)    :: A_inv1,A_inv2,A_inv3,A_inv4
        real, intent(in)    :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        end subroutine ib_solid_stress
        end interface

        interface
        subroutine ib_imBound( blockCount, blockList, timeEndAdv, dt)
        implicit none
        integer, INTENT(INOUT) :: blockCount
        integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList !blockCount
        real,    INTENT(IN) :: timeEndAdv,dt
        end subroutine
        end interface

        interface
        subroutine ib_ustar_solid(ustr, vstr, xms, Tau1,Tau2,Tau3,Tau4,&
                                  ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
        implicit none
        real, dimension(:,:,:), intent(inout) :: ustr, vstr
        real, dimension(:,:,:), intent(in)    :: xms,Tau1,Tau2,Tau3,Tau4
        real, intent(in)    :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        end subroutine ib_ustar_solid
        end interface

end module ib_viscoElastic_interface
