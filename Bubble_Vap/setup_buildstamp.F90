       subroutine setup_buildstamp (s_stamp_str, b_stamp_str, str_len)
       implicit none
       integer :: str_len
       character(len=str_len) :: s_stamp_str, b_stamp_str
       s_stamp_str = 'Fri Aug 27 06:08:20 2021'
       b_stamp_str = 'Fri Apr 1 10:31:56 2022'
       return
       end subroutine
      
       subroutine setup_systemInfo (system_str, str_len)
       integer :: str_len
       character(len=str_len) :: system_str
       system_str = 'Linux&
& pfe26&
& 4.12.14-122.80.1.20210720-nasa&
& #1 SMP Thu Jul 15 22:27:34 UTC 2021 (d956375)&
& x86_64'
       return
       end subroutine
      
