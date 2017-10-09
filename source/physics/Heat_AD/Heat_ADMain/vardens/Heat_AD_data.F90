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

     real, save    :: ht_Qmic

end module Heat_AD_data
