       subroutine setup_buildstamp (s_stamp_str, b_stamp_str, str_len)
       implicit none
       integer :: str_len
       character(len=str_len) :: s_stamp_str, b_stamp_str
       s_stamp_str = 'Thu Jun 10 12:32:34 2021'
       b_stamp_str = 'Thu Jun 10 12:41:00 2021'
       return
       end subroutine
      
       subroutine setup_systemInfo (system_str, str_len)
       integer :: str_len
       character(len=str_len) :: system_str
       system_str = 'Linux&
& pfe23&
& 4.12.14-122.66.2.20210415-nasa&
& #1 SMP Mon Apr 12 10:54:37 UTC 2021 (8fc5925)&
& x86_64'
       return
       end subroutine
      
