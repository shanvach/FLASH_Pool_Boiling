       subroutine setup_buildstamp (s_stamp_str, b_stamp_str, str_len)
       implicit none
       integer :: str_len
       character(len=str_len) :: s_stamp_str, b_stamp_str
       s_stamp_str = 'Sat Apr 1 21:53:29 2023'
       b_stamp_str = 'Sat Apr 1 22:00:14 2023'
       return
       end subroutine
      
       subroutine setup_systemInfo (system_str, str_len)
       integer :: str_len
       character(len=str_len) :: system_str
       system_str = 'Linux&
& pfe26&
& 3.10.0-1160.88.1.1chaos.ch6.x86_64&
& #1 SMP Fri Mar 17 14:57:36 PDT 2023&
& x86_64'
       return
       end subroutine
      
