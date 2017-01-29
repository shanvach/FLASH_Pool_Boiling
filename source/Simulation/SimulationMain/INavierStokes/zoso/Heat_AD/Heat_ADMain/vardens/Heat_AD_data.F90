module Heat_AD_data
     
     implicit none
     real, save :: ht_Pr_liq ! Prandtl Number for Liquid Phase
     real, save :: ht_Pr_gas ! Prandtl Number for Gas Phase
     real, save :: ht_Re_liq ! Re for Liquid Phase
     real, save :: ht_Re_gas ! Re for Gas Phase

     real, save :: ht_Nu ! Nusselt Number
     real, save :: ht_Bi ! Biot Number
     real, save :: ht_Twall_high, ht_Twall_low ! Non Dimensional wall Temp
     real, save :: ht_Tsat ! Saturation Temperature of the interface

end module Heat_AD_data
