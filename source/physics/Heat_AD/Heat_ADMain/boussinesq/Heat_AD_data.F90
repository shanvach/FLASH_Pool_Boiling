module Heat_AD_data
     
     implicit none
     real, save :: ht_invsqrtRaPr ! 1/SQRT(RaPr)
     real, save :: ht_source_term ! Internal Heat [0 || 1]
     real, save :: ht_Txl_value, ht_Txr_value, ht_Tyl_value, ht_Tyr_value
     real, save :: ht_Tzl_value, ht_Tzr_value
     integer, save :: ht_Txl_type, ht_Txr_type, ht_Tyl_type, ht_Tyr_type
     integer, save :: ht_Tzl_type, ht_Tzr_type

end module Heat_AD_data
