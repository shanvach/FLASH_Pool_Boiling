subroutine ImBound_setData()

!#ifdef INS_BOUSSINESQ
!  use IncompNS_data, only : ins_invsqrtRa_Pr
!#else
  use IncompNS_data, only : ins_invRe
!#endif

  use ImBound_data,  only : ib_nu 

  implicit none

!#ifdef INS_BOUSSINESQ
!  ib_nu = ins_invsqrtRa_Pr
!#else
  ib_nu = ins_invRe
!#endif

  return
end subroutine ImBound_setData
