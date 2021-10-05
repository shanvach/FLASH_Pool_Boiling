subroutine Heat_AD_finalize()

      use Heat_AD_data, only: ht_ibx, ht_iby, ht_ibz, ht_ibT, ht_ibNu

      implicit none

      deallocate(ht_ibx, ht_iby, ht_ibz, ht_ibT, ht_ibNu)

end subroutine Heat_AD_finalize
