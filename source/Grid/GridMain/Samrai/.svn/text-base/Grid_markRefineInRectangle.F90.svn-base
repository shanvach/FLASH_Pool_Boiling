!!****if* source/Grid/GridMain/Samrai/Grid_markRefineInRectangle
!!
!! NAME
!!  Grid_markRefineInRectangle
!!
!! SYNOPSIS
!!  Grid_markRefineInRectangle(real : ilb, irb, 
!!                        real : jlb, jrb, 
!!                        real : zlb, zrb, 
!!                        integer : lref, 
!!                        integer : contained)
!!
!! PURPOSE
!!  Refine blocks containing any points within a given rectangular 
!!  region having lower left coordinate (xlb,ylb,zlb) and upper right
!!  coordinate (xrb,yrb,zrb).  "Rectangular" is interpreted on a
!!  dimension-by-dimension basis: the region is an interval, rectangle,
!!  or rectangular parallelipiped in 1/2/3D Cartesian geometry; the
!!  rectangular cross-section of a rectangular torus in 2D axisymmetric
!!  (r-z) cylindrical geometry; an annular wedge in 2D polar (r-theta)
!!  cylindrical geometry; or an annulus in 1D spherical (r) geometry.
!!  Either blocks are brought up to a specific level of refinement or
!!  each block is refined once.

!!
!! ARGUMENTS
!!   ilb, irb,
!!   jlb, jrb, Bounding coords of the rectangular region to refine
!!   klb, krb
!!   lref      If > 0, bring all blocks to this level of refinement.
!!             If <= 0, refine all blocks once.
!!   contained If /= 0, refine only blocks completely contained within
!!             the rectangle; otherwise refine blocks with any overlap.
!!
!!***

subroutine Grid_markRefineInRectangle(ilb, irb, jlb, jrb, klb, krb, lref, contained)

!-------------------------------------------------------------------------------


implicit none
  real, intent(IN)    :: ilb, irb, jlb, jrb, klb, krb
  integer, intent(IN) :: lref, contained

  
  return
end subroutine Grid_markRefineInRectangle
