!!****ih* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhLocalInterface
!!
!! This is the header file for the Barnes-Hut tree solver that defines
!! additional interfaces private to the BHTree implementation.
!!***
Module gr_bhLocalInterface

  interface
    integer function gr_bhGetTreePos(level, mi)
      use gr_bhData, ONLY: gr_bhTreeLevels
      implicit none
      integer,intent(in) :: level, mi(1:gr_bhTreeLevels)
    end function gr_bhGetTreePos
  end interface

  interface gr_bhILContrib
    subroutine gr_bhILContrib(block, tr, cpu, temp, Phi, cellcnt, nodecnt)
      use gr_bhData, ONLY: gr_bhTreeBS
      implicit none
      integer,intent(in) :: block, tr, cpu, temp
      real,intent(INOUT),dimension(1:gr_bhTreeBS, 1:gr_bhTreeBS, 1:gr_bhTreeBS) :: Phi
      real,intent(out)   :: cellcnt, nodecnt
    end subroutine gr_bhILContrib
  end interface


  interface gr_bhLeafContrib
    subroutine gr_bhLeafContrib(block, tr, cpu, Phi, cellcnt, nodecnt)
      use gr_bhData, ONLY: gr_bhTreeBS
      implicit none
      integer,intent(in) :: block, tr, cpu
      real,intent(out)   :: cellcnt, nodecnt
      real, intent(INOUT), dimension(1:gr_bhTreeBS, 1:gr_bhTreeBS, 1:gr_bhTreeBS) :: Phi
    end subroutine gr_bhLeafContrib
  end interface

  interface gr_bhParentContrib
    subroutine gr_bhParentContrib(block, tr, cpu, Phi)
      use gr_bhData, ONLY: gr_bhTreeBS
      implicit none
      integer,intent(in) :: block, tr, cpu
      real, intent(INOUT), dimension(1:gr_bhTreeBS, 1:gr_bhTreeBS, 1:gr_bhTreeBS) :: Phi
    end subroutine gr_bhParentContrib
  end interface

end Module gr_bhLocalInterface


