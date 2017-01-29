

subroutine IncompNSstats()

  use Driver_data, only : dr_simTime
  use IncompNS_data, only : ins_meshMe,ins_nstep
  use IncompNSstats_data, only : instats_n,instats_IntervalStep,instats_start_time
  use instats_interface, only : instats_velp_timeavg, instats_Restresses_timeavg 

  use gr_interface, ONLY : gr_findMean

  implicit none
#include "constants.h"
#include "Flash.h"

  logical :: stats_flg

  real :: w,w2

  stats_flg = (MOD(ins_nstep,instats_IntervalStep) .eq. 0)

  if (stats_flg) then

    if (dr_simTime .ge. instats_start_time) then

    ! Velocities Time average
    call instats_velp_timeavg(instats_n)

    ! Reynolds Stresses Time Average
    call instats_Restresses_timeavg(instats_n)

    ! etc.

    ! Update ensemble number:
    instats_n = instats_n + 1

    ! Monitor w vel:
    call gr_findMean(DUST_VAR,2,.false.,w)  ! <w> vel that has been interpolated to cell center in instats_vep_timeavg
    call gr_findMean(WWAV_VAR,2,.false.,w2) ! <w*w> + <w>*<w>   

    if (ins_meshMe .eq. MASTER_PE) then
       write(*,*) "Time averaged statistics Computed, ensemble num=",instats_n
       write(*,*) "Time+domain averaged <w^2>=",w2-w*w,", where mean velocity <w>=",w
    endif

    else

      if (ins_meshMe .eq. MASTER_PE) then
        write(*,*) ' '
        write(*,*) &
        'Statistics delayed, simTime=',dr_simTime,' less than Stats start time=',instats_start_time
      endif

    endif

  endif

  return

end subroutine

