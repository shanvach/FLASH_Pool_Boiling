       subroutine setup_buildstamp (s_stamp_str, b_stamp_str, str_len)
       implicit none
       integer :: str_len
       character(len=str_len) :: s_stamp_str, b_stamp_str
       s_stamp_str = 'Sun Aug 13 20:07:36 2023'
       b_stamp_str = 'Wed Apr 24 15:53:34 2024'
       return
       end subroutine
      
       subroutine setup_systemInfo (system_str, str_len)
       integer :: str_len
       character(len=str_len) :: system_str
       system_str = 'Linux&
& pfe25&
& 4.18.0-477.15.1.1toss.t4.x86_64&
& #1 SMP Fri Jul 14 14:04:42 PDT 2023&
& x86_64'
       return
       end subroutine
      
