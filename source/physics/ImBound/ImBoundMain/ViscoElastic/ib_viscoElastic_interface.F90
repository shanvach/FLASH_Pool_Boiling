module ib_viscoElastic_interface

        !interface
        !subroutine ib_dynamic_grid_advection(sd,stest,&
        !                                     ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
        !implicit none
        !real, dimension(:,:,:), intent(inout) :: sn,adf,adfx,adfy
        !real, dimension(:,:,:), intent(in)    :: sd,stest
        !integer, dimension(:),  intent(in)    :: cpt
        !real, intent(in)    :: dx,dy,dz
        !integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        !end subroutine ib_dynamic_grid_advection
        !end interface

        interface
        subroutine ib_levelset_constantprojection(s,so,u,v,ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz) 
        implicit none
        real, dimension(:,:,:), intent(inout) :: s,so
        real, dimension(:,:,:), intent(in)    :: u,v
        real, intent(in)    :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        end subroutine ib_levelset_constantprojection
        end interface

        interface
        subroutine ib_levelset_linearprojection(s,so,sn,u,v,ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
        implicit none
        real, dimension(:,:,:), intent(inout) :: s,so
        real, dimension(:,:,:), intent(in)    :: u,v
        real, dimension(:,:,:), intent(in)    :: sn
        real, intent(in)    :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        end subroutine ib_levelset_linearprojection
        end interface

        interface
        subroutine ib_solid_stress(sd,sX,sY,Tau1,Tau2,Tau3,Tau4,&
                                   A1,A2,A3,A4,AT1,AT2,AT3,AT4,&
                                   A_inv1,A_inv2,A_inv3,A_inv4,&
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
        subroutine ib_advect( blockCount, blockList, timeEndAdv, dt)
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

        interface 
        subroutine ib_dynamic_grid_directional_derivative(sd,stest,adfx,adfy,sn,&
                                             ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
        implicit none
          real, dimension(:,:,:), intent(inout) :: sn,adfx,adfy
          real, dimension(:,:,:), intent(in)    :: sd,stest
          real, intent(in)    :: dx,dy,dz
          integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        end subroutine ib_dynamic_grid_directional_derivative
        end interface

        interface
        subroutine ib_dynamic_grid_normal_vector(sd,stest,adf,adfx,adfy,&
                                             ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
        implicit none
          real, dimension(:,:,:), intent(inout) :: adfx,adfy
          real, dimension(:,:,:), intent(in)    :: sd,stest
          real, dimension(:,:,:), intent(in)    :: adf 
          real, intent(in)    :: dx,dy,dz
          integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        end subroutine ib_dynamic_grid_normal_vector
        end interface

        interface
        subroutine ib_dynamic_grid_retain_inside(sd,sn0,sn,&
                                             ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
        implicit none
          real, dimension(:,:,:), intent(inout) :: sn0,sn
          real, dimension(:,:,:), intent(in)    :: sd
          real, intent(in)    :: dx,dy,dz
          integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        end subroutine
        end interface

        interface
        subroutine ib_redistance_PM(s,rho1, rho2, xmu1, xmu2, xmus,&
                                       rho, xmu, xms, blockID,     &
                                       phi, s1, dns,               &
                                       ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
          implicit none
            real, dimension(:,:,:), intent(inout) :: s, rho, xmu, xms
            real, intent(in)    :: rho1, rho2, xmu1, xmu2, xmus
            real, intent(in)    :: dx,dy,dz
            integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
            integer, intent(in) :: blockID
            real, dimension(:,:,:), intent(inout) :: phi, s1, dns
        end subroutine
        end interface

        interface
        subroutine ib_solid_interface_advection(sd,sX,sY,&
                                     ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
          implicit none
          real, dimension(:,:,:), intent(inout) :: sd,sX,sY
          real, intent(in)    :: dx,dy,dz
          integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        end subroutine
        end interface

end module ib_viscoElastic_interface
