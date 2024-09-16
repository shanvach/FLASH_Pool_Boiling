       subroutine setup_buildstamp (s_stamp_str, b_stamp_str, str_len)
       implicit none
       integer :: str_len
       character(len=str_len) :: s_stamp_str, b_stamp_str
       s_stamp_str = 'Thu Mar 23 14:10:15 2023'
       b_stamp_str = 'Wed Jul 26 13:50:08 2023'
       return
       end subroutine
      
       subroutine setup_systemInfo (system_str, str_len)
       integer :: str_len
       character(len=str_len) :: system_str
       system_str = 'Linux&
& pfe25&
& 3.10.0-1160.83.1.1chaos.ch6.x86_64&
& #1 SMP Tue Jan 24 17:36:15 PST 2023&
& x86_64'
       return
       end subroutine
      
