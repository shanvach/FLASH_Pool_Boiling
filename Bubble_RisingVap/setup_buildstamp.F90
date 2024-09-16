       subroutine setup_buildstamp (s_stamp_str, b_stamp_str, str_len)
       implicit none
       integer :: str_len
       character(len=str_len) :: s_stamp_str, b_stamp_str
       s_stamp_str = 'Fri Sep 3 13:17:45 2021'
       b_stamp_str = 'Mon Apr 18 19:02:59 2022'
       return
       end subroutine
      
       subroutine setup_systemInfo (system_str, str_len)
       integer :: str_len
       character(len=str_len) :: system_str
       system_str = 'Linux&
& pfe22&
& 4.12.14-122.83.1.20210810-nasa&
& #1 SMP Tue Aug 3 08:37:22 UTC 2021 (c86c48c)&
& x86_64'
       return
       end subroutine
      
