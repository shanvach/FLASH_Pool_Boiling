subroutine Heat_GFMstencil_o1(Tg,Ti,Tsat,th)

      implicit none

      real, intent(inout) :: Tg
      real, intent(in) :: Ti,th,Tsat

      Tg = (Tsat + (th-1)*Ti)/th


end subroutine Heat_GFMstencil_o1

subroutine Heat_GFMstencil_o2(Tg,Ti,Tip,Tsat,th)

      implicit none

      real, intent(inout) :: Tg
      real, intent(in) :: Ti,Tip,th,Tsat

      Tg = (2*Tsat + (2*th*th - 2)*Ti + (-th*th + th)*Tip)/(th + th*th)

end subroutine Heat_GFMstencil_o2
