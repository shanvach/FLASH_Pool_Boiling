!       Polynomial interpolation routine from Numerical Recipes


      subroutine polint(xa,ya,n,x,y,dy)
      
      implicit none

!..
!..given arrays xa and ya of length n and a value x, this routine returns a 
!..value y and an error estimate dy. if p(x) is the polynomial of degree n-1
!..such that ya = p(xa) ya then the returned value is y = p(x) 
!..
!..declare
      integer          n,nmax,ns,i,m
      parameter        (nmax=10)
      real             xa(n),ya(n),x,y,dy,c(nmax),d(nmax),dif,dift, & 
     &                 ho,hp,w,den
!..
!..find the index ns of the closest table entry; initialize the c and d tables
      ns  = 1
      dif = abs(x - xa(1))
      do 11 i=1,n
       dift = abs(x - xa(i))
       if (dift .lt. dif) then
        ns  = i
        dif = dift
       end if
       c(i)  = ya(i)
       d(i)  = ya(i)
11    continue
!..
!..first guess for y
      y = ya(ns)
!..
!..for each column of the table, loop over the c's and d's and update them
      ns = ns - 1
      do 13 m=1,n-1
       do 12 i=1,n-m
        ho   = xa(i) - x
        hp   = xa(i+m) - x
        w    = c(i+1) - d(i)
        den  = ho - hp
        if (den .eq. 0.0) stop ' 2 xa entries are the same in polint'
        den  = w/den
        d(i) = hp * den
        c(i) = ho * den
12     continue
!..
!..after each column is completed, decide which correction c or d, to add
!..to the accumulating value of y, that is, which path to take in the table
!..by forking up or down. ns is updated as we go to keep track of where we
!..are. the last dy added is the error indicator.
       if (2*ns .lt. n-m) then
        dy = c(ns+1)
       else
        dy = d(ns)
        ns = ns - 1
       end if
       y = y + dy
13    continue
      return
      end
