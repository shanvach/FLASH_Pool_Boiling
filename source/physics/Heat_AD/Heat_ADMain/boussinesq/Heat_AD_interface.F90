module Heat_AD_interface

    implicit none
        
    interface
       subroutine Heat_Solve3d(T_p, T_o, u, v, w, dt, dx, dy, dz, &
                               ix1, ix2, jy1, jy2, kz1, kz2)
         implicit none
         real, dimension(:,:,:), intent(inout) :: T_p
         real, dimension(:,:,:), intent(in) :: T_o
         real, dimension(:,:,:), intent(in) :: u,v,w
         real, intent(in) :: dt, dx, dy, dz 
         integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
       end subroutine Heat_Solve3d
    end interface

    interface
       subroutine Heat_Solve2d(T_p, T_o, u, v, dt, dx, dy, dz, ix1, ix2, jy1, jy2)
         implicit none
         real, dimension(:,:,:), intent(inout) :: T_p
         real, dimension(:,:,:), intent(in) :: T_o
         real, dimension(:,:,:), intent(in) :: u,v
         real, intent(in) :: dt, dx, dy, dz 
         integer, intent(in) :: ix1, ix2, jy1, jy2
       end subroutine Heat_Solve2d
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
