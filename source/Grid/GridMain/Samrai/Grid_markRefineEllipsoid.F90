!!****if* source/Grid/GridMain/Samrai/Grid_markRefineEllipsoid
!!  
!! NAME 
!!  Grid_markRefineEllipsoid 
!!  
!! SYNOPSIS 
!!  Grid_markRefineEllipsoid(real : x, 
!!                             real : y, 
!!                             real : z, 
!!                             real : a1, 
!!                             real : a2, 
!!                             real : a3, 
!!                             integer : lref) 
!!  
!! PURPOSE 
!!  Refine all blocks containing points on an ellipsoidal surface centered on
!!  (ic,jc,kc) with semimajor axes (a1,a2,a3).  Either blocks are brought up to
!!  a specific level of refinement or each block is refined once.  
!!  
!! ARGUMENTS 
!!  ic, jc, kc   Center of the ellipsoid
!!               (Coordinates for nonexistent dimensions are ignored.)
!!  a1, a2, a3   Semimajor axes of the ellipsoid
!!               (Axes for nonexistent dimensions are ignored.)
!!  lref         If > 0, bring all blocks to this level of refinement.
!!               If <= 0, refine all blocks once.
!!***  

subroutine Grid_markRefineEllipsoid(ic, jc, kc, a1, a2, a3, lref)


  implicit none

  real, intent(IN)      :: ic, jc, kc, a1, a2, a3
  integer, intent(IN)   :: lref

  
  return
end subroutine Grid_markRefineEllipsoid

