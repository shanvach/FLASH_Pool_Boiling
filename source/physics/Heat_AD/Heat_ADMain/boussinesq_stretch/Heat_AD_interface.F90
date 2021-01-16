module Heat_AD_interface

    implicit none
        
    interface
       subroutine ht_rhs3d(T, Trhs, u, v, w, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2)
         implicit none
         real, dimension(:,:,:), intent(in)  :: T
         real, dimension(:,:,:), intent(out) :: Trhs
         real, dimension(:,:,:), intent(in)  :: u, v, w
         real, dimension(:,:),   intent(in)  :: dx, dy, dz 
         integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
       end subroutine ht_rhs3d
    end interface

    interface
       subroutine ht_rhs2d(T, Trhs, u, v, dx, dy, ix1, ix2, jy1, jy2)
         implicit none
         real, dimension(:,:,:), intent(in)  :: T
         real, dimension(:,:,:), intent(out) :: Trhs
         real, dimension(:,:,:), intent(in)  :: u, v
         real, dimension(:,:),   intent(in)  :: dx, dy 
         integer, intent(in) :: ix1, ix2, jy1, jy2
       end subroutine ht_rhs2d
    end interface

    interface
       subroutine Heat_Solve(T, Trhs, Told, dt, ix1, ix2, jy1, jy2, kz1, kz2, gama, rho)
         implicit none
         real, dimension(:,:,:), intent(inout) :: T
         real, dimension(:,:,:), intent(in)    :: Trhs, Told
         real,    intent(in) :: dt, gama, rho
         integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
       end subroutine Heat_Solve
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

end module Heat_AD_interface
