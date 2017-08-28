module Heat_AD_interface

    implicit none
        
    interface
       subroutine Heat_Solve(T_p, T_o, T_rhs, dt, ix1, ix2, jy1, jy2, kz1, kz2, T_res)
         implicit none
         real, dimension(:,:,:), intent(inout) :: T_p
         real, dimension(:,:,:), intent(in) :: T_o
         real, dimension(:,:,:), intent(in) :: T_rhs
         real, intent(in) :: dt
         integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
         real, intent(inout) :: T_res
       end subroutine Heat_Solve
    end interface

    interface
       subroutine Heat_RHS_upwind(T_rhs, T_o, u, v, dx, dy, dz, inRe, ix1, ix2, jy1, jy2,&
                           rho1x,rho2x,rho1y,rho2y,alph,pf,s,mdot,nrmx,nrmy,smrh,curv)
         implicit none
         real, dimension(:,:,:), intent(inout) :: T_rhs
         real, dimension(:,:,:), intent(in) :: T_o
         real, dimension(:,:,:), intent(in) :: u,v
         real, intent(in) :: dx, dy, dz, inRe
         integer, intent(in) :: ix1, ix2, jy1, jy2
         real, dimension(:,:,:),intent(in) :: rho1x,rho2x,rho1y,rho2y,alph
         real, dimension(:,:,:),intent(in) :: pf,s,mdot,nrmx,nrmy,smrh,curv
       end subroutine Heat_RHS_upwind
    end interface

    interface
       subroutine Heat_RHS_central(T_rhs, T_o, u, v, dx, dy, dz, inRe, ix1, ix2, jy1, jy2,&
                           rho1x,rho2x,rho1y,rho2y,alph,pf,s,mdot,nrmx,nrmy,smrh,curv)
         implicit none
         real, dimension(:,:,:), intent(inout) :: T_rhs
         real, dimension(:,:,:), intent(in) :: T_o
         real, dimension(:,:,:), intent(in) :: u,v
         real, intent(in) :: dx, dy, dz, inRe
         integer, intent(in) :: ix1, ix2, jy1, jy2
         real, dimension(:,:,:),intent(in) :: rho1x,rho2x,rho1y,rho2y,alph
         real, dimension(:,:,:),intent(in) :: pf,s,mdot,nrmx,nrmy,smrh,curv
       end subroutine Heat_RHS_central
    end interface

    interface
       subroutine Heat_RHS_weno3(T_rhs, T_o, u, v, dx, dy, dz, inRe, ix1, ix2, jy1, jy2,&
                           rho1x,rho2x,rho1y,rho2y,alph,pf,s,mdot,nrmx,nrmy,smrh,curv)
         implicit none
         real, dimension(:,:,:), intent(inout) :: T_rhs
         real, dimension(:,:,:), intent(in) :: T_o
         real, dimension(:,:,:), intent(in) :: u,v
         real, intent(in) :: dx, dy, dz, inRe
         integer, intent(in) :: ix1, ix2, jy1, jy2
         real, dimension(:,:,:),intent(in) :: rho1x,rho2x,rho1y,rho2y,alph
         real, dimension(:,:,:),intent(in) :: pf,s,mdot,nrmx,nrmy,smrh,curv
       end subroutine Heat_RHS_weno3
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
      subroutine Heat_calGradT_central(Tnl,Tnv,T,s,pf,dx,dy,dz,ix1,ix2,jy1,jy2,nx,ny,mflg)
        implicit none
        real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
        real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny,mflg
        real, intent(in) :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2
      end subroutine Heat_calGradT_central
   end interface

   interface 
      subroutine Heat_calGradT(Tnl,Tnv,T,s,pf,dx,dy,dz,ix1,ix2,jy1,jy2,nx,ny,mflg)
        implicit none
        real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
        real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny,mflg
        real, intent(in) :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2
      end subroutine Heat_calGradT
   end interface

   interface
      subroutine Heat_extrapGradT(Tnl,Tnv,T,s,pf,dx,dy,dz,nx,ny,ix1,ix2,jy1,jy2,Tnl_res,Tnv_res,mflg)
        implicit none
        real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
        real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny,mflg
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
       subroutine Heat_RHS_3D(T_rhs, T_o, u, v, w, dx, dy, dz, inRe, ix1, ix2, jy1, jy2, kz1, kz2,&
                              rho1x,rho2x,rho1y,rho2y,rho1z,rho2z,alph,pf,s,mdot,nrmx,nrmy,nrmz,smrh,curv)
         implicit none
         real, dimension(:,:,:), intent(inout) :: T_rhs
         real, dimension(:,:,:), intent(in) :: T_o
         real, dimension(:,:,:), intent(in) :: u,v,w
         real, intent(in) :: dx, dy, dz, inRe
         integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
         real, dimension(:,:,:),intent(in) :: rho1x,rho2x,rho1y,rho2y,alph,rho1z,rho2z
         real, dimension(:,:,:),intent(in) :: pf,s,mdot,nrmx,nrmy,nrmz,smrh,curv
       end subroutine Heat_RHS_3D
   end interface

   interface
       subroutine Heat_RHS_3D_weno3(T_rhs, T_o, u, v, w, dx, dy, dz, inRe, ix1, ix2, jy1, jy2, kz1, kz2,&
                              rho1x,rho2x,rho1y,rho2y,rho1z,rho2z,alph,pf,s,mdot,nrmx,nrmy,nrmz,smrh,curv)
         implicit none
         real, dimension(:,:,:), intent(inout) :: T_rhs
         real, dimension(:,:,:), intent(in) :: T_o
         real, dimension(:,:,:), intent(in) :: u,v,w
         real, intent(in) :: dx, dy, dz, inRe
         integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
         real, dimension(:,:,:),intent(in) :: rho1x,rho2x,rho1y,rho2y,alph,rho1z,rho2z
         real, dimension(:,:,:),intent(in) :: pf,s,mdot,nrmx,nrmy,nrmz,smrh,curv
       end subroutine Heat_RHS_3D_weno3
   end interface

   interface 
      subroutine Heat_calGradT_3D(Tnl,Tnv,T,s,pf,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,nx,ny,nz,mflg)
        implicit none
        real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
        real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny,nz,mflg
        real, intent(in) :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
      end subroutine Heat_calGradT_3D
   end interface

   interface 
      subroutine Heat_calGradT_3D_central(Tnl,Tnv,T,s,pf,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,nx,ny,nz,mflg)
        implicit none
        real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
        real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny,nz,mflg
        real, intent(in) :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
      end subroutine Heat_calGradT_3D_central
   end interface

   interface
      subroutine Heat_extrapGradT_3D(Tnl,Tnv,T,s,pf,dx,dy,dz,nx,ny,nz,ix1,ix2,jy1,jy2,kz1,kz2,Tnl_res,Tnv_res,mflg)
        implicit none
        real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
        real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny,nz,mflg
        real, intent(in) :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        real, intent(out) :: Tnl_res,Tnv_res
      end subroutine Heat_extrapGradT_3D
   end interface

   interface
      subroutine Heat_GFMstencil_o1(Tg,Ti,Tsat,th)
        implicit none
        real, intent(inout) :: Tg
        real, intent(in) :: Ti,th,Tsat
      end subroutine

      subroutine Heat_GFMstencil_o2(Tg,Ti,Tip,Tsat,th)
        implicit none
        real, intent(inout) :: Tg
        real, intent(in) :: Ti,Tip,th,Tsat
      end subroutine
   end interface

end module Heat_AD_interface
