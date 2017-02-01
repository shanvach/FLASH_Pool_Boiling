module Heat_AD_interface

    implicit none
        
    interface
       subroutine Heat_Solve(T_p, T_o, T_rhs, dt, ix1, ix2, jy1, jy2,T_res)
         implicit none
         real, dimension(:,:,:), intent(inout) :: T_p
         real, dimension(:,:,:), intent(in) :: T_o
         real, dimension(:,:,:), intent(in) :: T_rhs
         real, intent(in) :: dt
         integer, intent(in) :: ix1, ix2, jy1, jy2
         real, intent(inout) :: T_res
       end subroutine Heat_Solve
    end interface

    interface
       subroutine Heat_RHS(T_rhs, T_o, u, v, dx, dy, dz, inRe, ix1, ix2, jy1, jy2,rho1x,rho2x,rho1y,rho2y,thco,cp,pf,s,mdot,nrmx,nrmy,smrh)
         implicit none
         real, dimension(:,:,:), intent(inout) :: T_rhs
         real, dimension(:,:,:), intent(in) :: T_o
         real, dimension(:,:,:), intent(in) :: u,v
         real, intent(in) :: dx, dy, dz, inRe
         integer, intent(in) :: ix1, ix2, jy1, jy2
         real, dimension(:,:,:),intent(in) :: rho1x,rho2x,rho1y,rho2y,thco,cp
         real, dimension(:,:,:),intent(in) :: pf,s,mdot,nrmx,nrmy,smrh
       end subroutine Heat_RHS
    end interface

    interface
      subroutine Heat_AD_init(blockCount,blockList)
      implicit none
      integer, INTENT(INOUT) :: blockCount
      integer, INTENT(INOUT) :: blockList(MAXBLOCKS)
      end subroutine Heat_AD_init
    end interface

    interface
        subroutine Heat_AD_finalize()
          implicit none
        end subroutine
    end interface

    interface
       subroutine Heat_AD(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)
         implicit none
         integer, INTENT(INOUT) :: blockCount
         integer, INTENT(IN) :: sweepOrder
         integer, INTENT(INOUT) :: blockList(MAXBLOCKS)
         real,    INTENT(IN) :: timeEndAdv, dt, dtOld
       end subroutine Heat_AD
    end interface 

   interface 
      subroutine Heat_calGradT(Tnl,Tnv,T,s,pf,dx,dy,dz,ix1,ix2,jy1,jy2,nx,ny)
        implicit none
        real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
        real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny
        real, intent(in) :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2
      end subroutine Heat_calGradT
   end interface

   interface
      subroutine Heat_extrapGradT(Tnl,Tnv,T,s,pf,dx,dy,dz,nx,ny,ix1,ix2,jy1,jy2,Tnl_res,Tnv_res)
        implicit none
        real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
        real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny
        real, intent(in) :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2
        real, intent(out) :: Tnl_res,Tnv_res
      end subroutine Heat_extrapGradT
   end interface

   interface
      subroutine Heat_calMdot(mdot,Tnl,Tnv,alpha_l,alpha_v,nx,ny,ix1,ix2,jy1,jy2,kz1,kz2)
        implicit none
        real, dimension(:,:,:), intent(inout) :: mdot
        real, dimension(:,:,:), intent(in) :: Tnl,Tnv,nx,ny
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        real, intent(in) :: alpha_l,alpha_v
      end subroutine Heat_calMdot
   end interface

   interface
       subroutine Heat_Solve_3D(T_p, T_o, T_rhs, dt, ix1, ix2, jy1, jy2,T_res)
         implicit none
         real, dimension(:,:,:), intent(inout) :: T_p
         real, dimension(:,:,:), intent(in) :: T_o
         real, dimension(:,:,:), intent(in) :: T_rhs
         real, intent(in) :: dt
         integer, intent(in) :: ix1, ix2, jy1, jy2
         real, intent(inout) :: T_res
       end subroutine Heat_Solve_3D
   end interface

   interface
       subroutine Heat_RHS_3D(T_rhs, T_o, u, v, dx, dy, dz, inRe, ix1, ix2, jy1, jy2,rho1x,rho2x,rho1y,rho2y,thco,cp,pf,s,mdot,nrmx,nrmy,smrh)
         implicit none
         real, dimension(:,:,:), intent(inout) :: T_rhs
         real, dimension(:,:,:), intent(in) :: T_o
         real, dimension(:,:,:), intent(in) :: u,v
         real, intent(in) :: dx, dy, dz, inRe
         integer, intent(in) :: ix1, ix2, jy1, jy2
         real, dimension(:,:,:),intent(in) :: rho1x,rho2x,rho1y,rho2y,thco,cp
         real, dimension(:,:,:),intent(in) :: pf,s,mdot,nrmx,nrmy,smrh
       end subroutine Heat_RHS_3D
   end interface

   interface 
      subroutine Heat_calGradT_3D(Tnl,Tnv,T,s,pf,dx,dy,dz,ix1,ix2,jy1,jy2,nx,ny)
        implicit none
        real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
        real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny
        real, intent(in) :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2
      end subroutine Heat_calGradT_3D
   end interface

   interface
      subroutine Heat_extrapGradT_3D(Tnl,Tnv,T,s,pf,dx,dy,dz,nx,ny,ix1,ix2,jy1,jy2,Tnl_res,Tnv_res)
        implicit none
        real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
        real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny
        real, intent(in) :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2
        real, intent(out) :: Tnl_res,Tnv_res
      end subroutine Heat_extrapGradT_3D
   end interface

end module Heat_AD_interface
