subroutine  MACHINEEPSILON
  use Simulation_data, ONLY: eps
  implicit none 
  double percision :: MACHEPS=1.D0

  MACHEPS = MACHEPS / 2.D0
  !MACHEPS = 1.D0
  do while (1.D0 + MACHEPS / 2.D0 .NE. 1.D0) then
     MACHEPS = MACHEPS / 2.D0
     !  IF ( 1.D0 + MACHEPS / 2.D0 .EQ. 1.D0 ) GOTO 110
  end do
  !GO TO 100
  !110  CONTINUE
  eps= MACHEPS;
  print*, MACHEPS
  return
end subroutine MACHINEEPSILON
