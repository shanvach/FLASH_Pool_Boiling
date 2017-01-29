!       Binary search routine from Numerical Recipes


      subroutine hunt(xx,n,x,jlo) 

      implicit none

!.. 
!..given an array xx of length n and a value x, this routine returns a value 
!..jlo such that x is between xx(jlo) and xx(jlo+1). the array xx must be  
!..monotonic. j=0 or j=n indicates that x is out of range. on input jlo is a  
!..guess for the table entry, and should be used as input for the next hunt.  
!..one assumes the next table entry is relatively close to the old one. 
!.. 
!..declare 
      logical          ascnd 
      integer          n,jlo,jhi,inc,jm 
      real             xx(n),x 
!.. 
!..logical function definition 
!..will be true if we are ascending the order of the table; false if descending 
      ascnd = xx(n) .ge. xx(1) 
!.. 
!..horrid initial guess; go to bisection immediatly 
      if (jlo.le.0 .or. jlo.gt.n) then 
       jlo = 0 
       jhi = n+1 
       go to 3 
      end if 
!.. 
!..set the initial hunting increment 
      inc = 1 
!.. 
!..this is the hunt up section 
      if (x.ge.xx(jlo) .eqv. ascnd) then 
1      jhi = jlo + inc 
!.. 
!..if we are the end of the table then we are done hunting 
       if (jhi .gt. n) then 
        jhi = n + 1 
!.. 
!..otherwise replace the lower limit, double the increment and try again 
       else if (x.ge.xx(jhi) .eqv. ascnd) then 
        jlo = jhi 
        inc = inc + inc 
        go to 1 
       end if 
!.. 
!..this is the hunt down section 
      else 
       jhi = jlo 
2      jlo = jhi - inc 
!.. 
!..if we are the end of the table then we are done hunting 
       if (jlo.lt.1) then 
        jlo = 0 
!.. 
!..otherwise replace the upper limit, double the increment and try again 
       else if (x.lt.xx(jlo) .eqv. ascnd) then 
        jhi = jlo 
        inc = inc + inc 
        go to 2 
       end if 
      end if 
!.. 
!..this is the final bisection phase 
3     if (jhi-jlo .eq. 1) then 
       if (x .eq. xx(n)) jlo = n - 1 
       if (x .eq. xx(1)) jlo = 1 
       return 
      end if 
      jm = (jhi+jlo)/2 
      if (x .ge. xx(jm) .eqv. ascnd) then 
       jlo = jm 
      else 
       jhi = jm 
      end if 
      go to 3 
      end 
