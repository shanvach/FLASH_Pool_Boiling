module Heat_AD_data
     
     implicit none
     real, save :: ht_Pr ! Prandtl Number for Liquid Phase

     real, save :: ht_Nu ! Nusselt Number
     real, save :: ht_Bi ! Biot Number
     real, save :: ht_Twall_high, ht_Twall_low ! Non Dimensional wall Temp
     real, save :: ht_Tsat ! Saturation Temperature of the interface

     real, save :: ht_coeff(3,3)
     real, save :: ht_L
     real, save :: ht_St

     integer, save :: ht_hfit

     real, save    :: ht_AMR_specs(2)

     real, save    :: ht_qmic, ht_fmic

     real, save    :: ht_dxmin

     real, save    :: ht_Ab, ht_Bb, ht_Cb

     real, save    :: ht_tWait

     real, save    :: ht_psi

     real, save    :: ht_Nu_l, ht_Nu_t

     real, save    :: ht_Ra = 0

     real, save    :: ht_Tnuc

     real, save    :: ht_convel

     logical, save :: ht_microFlg

     real, allocatable, save :: ht_ibx(:), ht_iby(:), ht_ibz(:), ht_ibT(:), ht_ibNu(:)

     real :: ht_ibNu_t, ht_ibNu_l

     integer, save :: ht_hflux_counter, ht_ibhflux_counter

     logical, save :: ht_hflux_flag

end module Heat_AD_data
