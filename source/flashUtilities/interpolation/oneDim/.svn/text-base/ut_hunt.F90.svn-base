!!****if* source/flashUtilities/interpolation/oneDim/ut_hunt
!!
!!  NAME
!!    ut_hunt  
!!  SYNOPSIS 
!!    call ut_hunt(xx, n, x, low)
!!  FUNCTION 
!!    Given an array xx of length n and a value x, this routine returns a value 
!!    low such that x is between xx(low) and xx(low+1). the array xx must be  
!!    monotonic. j=0 or j=n indicates that x is out of range. on input low is a  
!!    guess for the table entry, and should be used as input for the next hunt.  
!!  INPUTS
!!      xx - real, sorted, array of dimension n to search through
!!      n  - size of xx
!!      x  - value to search for
!!      low - largest value in [1,n] s.t. xx[low] < x
!!  RESULT
!!      most errors return low = 1
!!
!!  HISTORY
!!
!!    This was in FLASH2 as setups/nova/non-numrecipies/hunt.F90.
!!    Changed name to ut_hunt         - KW 2007-06-07
!!    Implemented documented use of 'low' value passed in - KW 2007-06-07
!!***
      subroutine ut_hunt(xx,n,x,low) 

      implicit none

      integer, INTENT(in)            :: n
      real, INTENT(in), DIMENSION(n) :: xx
      real, INTENT(in)               :: x
      integer, INTENT(inout)         :: low

! dirn < 0 for decreasing, >0 for increasing

      real                           :: dirn
      integer                        :: high, middle

      dirn = xx(n)-xx(1)

! this should never, ever happen.
      if (dirn == 0.) then
         print *,'ut_hunt: array xx() is constant!'
         low = 1
         return
      endif

! same here
      if (n < 2) then
         print *,'list too small!'
         low = 1
         return
      endif

      if (dirn > 0) then
         if (low .ge. 1 .AND. low < n) then
            if (xx(low) .le. x .AND. x < xx(low+1)) then
               return
            end if
         end if
         low = 1 ; high = n
      else
         if (low .ge. 1 .AND. low < n) then
            if (xx(low+1) .lt. x .AND. x < xx(low)) then
               return
            end if
         end if
         low = n ; high = 1
      endif

! same here
      if (xx(low) > x ) then
         if (dirn > 0) low = 0
         return
      endif

      if (xx(high) < x) then
         if (dirn < 0) then
            low = 0
         else
            low = n
         end if
         return
      endif

      if (dirn > 0) then
searchinc:  do while (high-low > 1) 
          middle = (low+high)/2
          if (xx(middle) .LE. x) then
                        low = middle
              else
            high = middle
          endif

      enddo searchinc
      else
          high = 1
searchdec:  do while (low - high > 1) 
          middle = (low+high)/2
          if (xx(middle) .LE. x) then
                        low = middle
              else
            high = middle
          endif
      enddo searchdec
      endif

      if (xx(low) > x .or. x > xx(high)) then
         print *,'Bug alert! bug alert!'
      endif

      end subroutine ut_hunt
