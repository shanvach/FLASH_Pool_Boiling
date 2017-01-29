!!****ih* source/Grid/localAPI/gr_bhInterface
!!
!! This is the header file for the Barnes-Hut tree solver that defines its
!! interfaces private to the Grid unit.
!!***
Module gr_bhInterface

  interface gr_bhInit
    subroutine gr_bhInit()
    end subroutine gr_bhInit
  end interface

  interface gr_bhFinalize
    subroutine gr_bhFinalize()
    end subroutine gr_bhFinalize
  end interface

  interface gr_bhBlockRelationship
    integer function gr_bhBlockRelationship(block, tr, cpu)
      integer, intent(in) :: block, tr, cpu
    end function gr_bhBlockRelationship
  end interface

  interface gr_bhBuildTree
    subroutine gr_bhBuildTree(idensvar)
      integer, intent(in) :: idensvar
    end subroutine gr_bhBuildTree
  end interface

  interface gr_bhBuildTreeBlock
    subroutine gr_bhBuildTreeBlock(idensvar, block)
      integer,intent(in) :: idensvar
      integer,intent(in) :: block
    end subroutine gr_bhBuildTreeBlock
  end interface

  interface gr_bhComBlkProperties
    subroutine gr_bhComBlkProperties()
    end subroutine gr_bhComBlkProperties
  end interface

  interface gr_bhComParentTree
    subroutine gr_bhComParentTree(level)
      integer, intent(IN) :: level
    end subroutine gr_bhComParentTree
  end interface

  interface gr_bhDestroyTree
    subroutine gr_bhDestroyTree()
    end subroutine gr_bhDestroyTree
  end interface

  interface gr_bhEwald
    real function gr_bhEwald(x, y, z)
      real,intent(in)    :: x, y, z
    end function gr_bhEwald
  end interface

  interface
    subroutine gr_bhEwaldField()
    end subroutine gr_bhEwaldField
  end interface

  interface gr_bhExchangeTrees
    subroutine gr_bhExchangeTrees()
    end subroutine gr_bhExchangeTrees
  end interface

  interface gr_bhFindNeighbours
    subroutine gr_bhFindNeighbours()
    end subroutine gr_bhFindNeighbours
  end interface


  interface gr_bhGetTreeSize
    integer function gr_bhGetTreeSize(level)
      integer,intent(in) :: level
    end function gr_bhGetTreeSize
  end interface

  interface gr_bhInitTemplates
    subroutine gr_bhInitTemplates(write_arrs, max_cc_dist, max_cc_ind, max_cn_ind)
      integer , intent(in)  :: write_arrs
      integer , intent(inout) :: max_cc_dist, max_cc_ind, max_cn_ind
    end subroutine gr_bhInitTemplates
  end interface


  interface gr_bhPotential
    subroutine gr_bhPotential(idensvar, ipotvar)
      integer, intent(in) :: ipotvar, idensvar
    end subroutine gr_bhPotential
  end interface

  interface gr_bhPotentialBlock
    subroutine gr_bhPotentialBlock(block, idensvar, ipotvar)
      integer, intent(IN) :: block, idensvar, ipotvar
    end subroutine gr_bhPotentialBlock
  end interface


end Module gr_bhInterface


