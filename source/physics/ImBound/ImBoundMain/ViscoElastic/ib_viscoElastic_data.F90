!! source/physics/ImBound/ViscoElastic/ib_viscoElastic_data
!!
!! NAME
!!
!!  ib_viscoElastic_data
!!
!!
!! SYNOPSIS
!!
!!  MODULE ib_viscoElastic_data()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  This stores data and limiter functions that are specific to the Visco Elastic sub unit
!!
!!***
 
 
module ib_viscoElastic_data

  real, save :: ib_vis_rho1, ib_vis_rho2  
  real, save :: ib_vis_xmu1, ib_vis_xmu2, ib_vis_xmus 
  !ib_vis_rho1 = 1.d0
  !ib_vis_rho2 = 1.d0
  !ib_vis_xmu1 = 1.d0
  !ib_vis_xmu2 = 1.d0
  !ib_vis_xmus = 1.d0

end module ib_viscoElastic_data
