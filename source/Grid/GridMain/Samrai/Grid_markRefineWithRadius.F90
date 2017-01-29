!!****if* source/Grid/GridMain/Samrai/Grid_markRefineWithRadius
!!  
!! NAME 
!!  Grid_markRefineWithRadius 
!!  
!! SYNOPSIS 
!!  Grid_markRefineWithRadius(real : ic, 
!!                       real : jc, 
!!                       real : kc, 
!!                       real : radius, 
!!                       integer : lref) 
!!  
!! PURPOSE 
!!  Refine all blocks containing points at a given distance from a given point
!!  (ic,jc,kc).  Either blocks are brought up to a specific level of refinement
!!  or each block is refined once.  
!!  
!! ARGUMENTS 
!!  ic, jc, kc   Center of the interval/circle/sphere
!!               (Coordinates for nonexistent dimensions are ignored.)
!!  radius       Radius of the region 
!!  lref         If > 0, bring all blocks to this level of refinement.
!!               If <= 0, refine all blocks once.
!!  
!!  
!!***

subroutine Grid_markRefineWithRadius(ic, jc, kc, radius, lref)


implicit none
  real, intent(IN)      :: ic, jc, kc, radius
  integer, intent(IN)   :: lref

  
  return
end subroutine Grid_markRefineWithRadius
