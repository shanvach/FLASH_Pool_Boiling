c///////////////////////////////////////////////////////////////////////////////
c//
c// Program to compute the analytic solution to the Taylor-Sedov blastwave
c// problem. The initial density and energy are set equal to 1. The Eulerian
c// solution is then used to compute the Lagrangian solution of a tracer
c// particle within the flow, beginning at a radius r_particle.
c//
c// Ratio of specific heats gamma is set to 7/5 by default below.
c//
c// Written by Robert Fisher, 12/02/05.
c// Based on code written by Robert Fisher and Jeffrey Jauregui.
c//
c///////////////////////////////////////////////////////////////////////////////


#define NMAX 10000
#define NSTEPS 100
#define RMAX 1.0
#define TMAX 1.0


      program main

      implicit none

      integer n, i
      real*8 pi, t, gamma, r, dr, dt, lagrangian
      real*8 density (NMAX), velocity (NMAX), pressure (NMAX),
     &        mass (0:NMAX), mass_r, r_particle, r_final

c... Constants

      pi = 4. * atan (1.)
      gamma = 1.4

c... Solution grid spacing

      dr = RMAX / NMAX
      dt = TMAX / NSTEPS

c... Initial particle location

      r_particle = 0.2175

c... Loop over times

      do n = 1, NSTEPS

        t = n * dt        

        call taylorsedov_solution(t, gamma, density, velocity, pressure)

c... Integrate density distribution over radius to obtain cumulative
c...   mass profile of TS solution at time t

        mass (0) = 0.

        do i = 1, NMAX
          r        = i * dr
          mass (i) = mass (i-1) + 4. * pi * r**2. * dr * density (i)
        end do

c... Compute Lagrangian mass coordinate of particle based on its initial
c...  position r = particle at t = 0.

        lagrangian = 4./3. * pi * r_particle**3.  ! initial density = 1

c... Solve for final particle radius by locating the Lagrangian mass coordinate
c...   of the particle in the final Eulerian density distribution.

        i = 0
        mass_r = mass (i)
     
        do i = 1, NMAX
          mass_r = mass (i)
          if (mass_r .gt. lagrangian) then
            goto 10
          endif
        end do

 10    r_final = i * dr   

c... Output solution table (t, r).

       print *, t, ' ', r_final

c... table of density, velocity, pressure at solution time t

c       do i = 1, NMAX
c         print *, i * dr, ' ', density (i), ' ', velocity (i), ' ',
c     &               pressure (i)
c       end do

       end do ! loop over NSTEPS

       end

c..........................................................................
c...
c... Subroutine to compute analytic functions. Currently we use this
c...  to determine the solution to the Taylor-Sedov problem, using
c...  an iterative numerical method to determine Landau & Lifschtiz's 
c...  "critical" velocity.
c...
c... The method used is a simple bisection method (rtbis), which is
c...  always guaranteed to find a solution.
c...
c... The region near the origin is treated specially, because the 
c...  the TS problem is singular there.
c... 
c...........................................................................

      subroutine taylorsedov_solution (t, gamma, rho_out,
     &                                 v_out, p_out)

      implicit none

      integer n 
      real*8 t, r, rho_out (NMAX), v_out (NMAX), p_out (NMAX)
      real*8 density, pressure, velocity, rtbis, gamma
      real*8 rtbis_TS, func_TS

      real*8 rho1, rho2, u1, v2, p1, p2, E0, cap_R, beta
      real*8 cap_G, nu1, nu2, nu3, nu4, nu5
      real*8 fudge_V, fudge_rho, v_crit
      real*8 cap_V, cap_Z

c     WARNING: the analytic solution will only
c     be computed correctly for gamma .ge. 1.3
 
      if (gamma .lt. 1.3) then
        print *, "Taylor-Sedov solution not accurate for gamma < 1.3"
        stop
      endif

c     See Landau & Lifshitz "Fluid Mechanics",
c     2nd edition, pg. 403-406 for variable definitions

      rho1 = 1.D0
      p1 = 1.D0/(gamma-1.D0) * 1.D-8
      E0 = 1.D0

c     interpolate the constant beta in terms of known 
c     values, via a 2nd order polynomial fit.  At 
c     gamma = 1.2, 1.4, 5/3, beta values are 
c     0.89, 1.033, and 1.17, respectively

      beta = -4.3125D-1*gamma**2.D0 + 1.83624D0*gamma
     &         - 6.9250D-1

c     these are fudge factors for the velocity and 
c     density that ensure the functions are continuous 
c     when we go from the exact to approximate solution 
c     near the origin.  fudge_rho is determined from 
c     a 2nd order polynomial fit from values at gamma 
c     equal to 1.3, 1.4, and 1.6

      fudge_V = 5.D-1*(1.D0 + 1.D0/gamma) 
      fudge_rho = 9.4867D-02*gamma**2.D0
     &   -1.62340D-1*gamma + 3.52397D-1


c     some constants

      nu1 = -1.D0*(13.D0*gamma**2.D0 - 7.D0*gamma + 12.D0)
     &        /((3.D0*gamma-1.D0)*(2.D0*gamma+1.D0))
      nu2 = 5.D0*(gamma-1.D0)/(2.D0*gamma+1.D0)
      nu3 = 3.D0/(2.D0 * gamma + 1.D0)
      nu4 = (13.D0 * gamma**2.D0 - 7.D0 * gamma + 12.D0)/
     &       ((3.D0*gamma-1.D0)*(2.D0*gamma + 1.D0)*(2.D0-gamma))
      nu5 = 2.D0 / (gamma - 2.D0)

c     cap_R is the position of the shock front at time t, 
c     u1 is its propagation speed

      cap_R = beta*((E0 * t**2.D0)/rho1)**2.D-1
      u1 = 2.D0 * beta * E0**2.D-1 / (5.D0 * rho1**2.D-1 * t**6.D-1)

      v2 = 2.D0 * u1 / (gamma+1)
      rho2 = rho1 * (gamma+1.D0) / (gamma-1.D0)
      p2 = 2.D0 * rho1 * u1**2.D0 / (gamma+1.D0)

      do n = 1, NMAX

              r = (RMAX / NMAX) * n 

c      if r is outside the shock front, do not solve 
c      for anything

              if (r .gt. cap_R) then
                density = rho1
                pressure = p1
                velocity = 0.D0
c             inside the shock front:
              else

c             if the self-similarity variable is greater than 0.2, 
c             just compute the solution normally

              if(r / cap_R .gt. 2.D-1) then
                cap_V = rtbis_TS(1.D0/gamma + 1.D-10,
     &                  5.D0/(3.D0*gamma-1.D0)-1.D-2,
     &                  r/cap_R, gamma, nu1, nu2, 1.D-16)

                cap_Z = gamma * (gamma - 1.D0)*(1.D0 - cap_V)
     &                   * cap_V**2.D0 / (2.D0 * (gamma * cap_V - 1.D0))

                cap_G = (gamma+1)/(gamma-1)
     &               *((gamma+1)/(gamma-1)*(gamma*cap_V-1))**nu3
     &               * ((gamma+1)/(7-gamma)*(5-(3*gamma-1)*cap_V))**nu4
     &               * ((gamma+1)/(gamma-1)*(1-cap_V))**nu5

                pressure = (4.D0 * r**2.D0 
     &              * cap_Z * rho1*cap_G)
     &              / (25.D0 * t**2.D0 * gamma)

c             inside r/cap_R = 0.2:
              else
c             use the approximate functional form of cap_V
                  cap_V = fudge_V*2.5D0 * t * v2 / cap_R

c                 determine the "critical" value of cap_V:
c                 That is, its value at r/cap_R = 0.2

                  v_crit = rtbis_TS(1.D0/gamma + 1.D-10,
     &                      5.D0/(3.D0*gamma-1.D0)-1.D-2,
     &                      2.D-1, gamma, nu1, nu2, 1.D-16)

c                 compute cap_Z and cap_G using v_crit

                  cap_Z = gamma * (gamma - 1.D0)*(1.D0 - v_crit)
     &                    * v_crit**2.D0 /
     &                       (2.D0 * (gamma * v_crit - 1.D0))
                  cap_G = (gamma+1)/(gamma-1)
     &              * ((gamma+1)/(gamma-1)*(gamma*v_crit-1))**nu3
     &              * ((gamma+1)/(7-gamma)*(5-(3*gamma-1)*v_crit))**nu4
     &              * ((gamma+1)/(gamma-1)*(1-v_crit))**nu5

c                 everywhere inside r/cap_R =0.2, set the pressure to 
c                 its value at 0.2
                  pressure = 4.D0*(2.D-1*cap_R)**2.D0*cap_Z*rho1*cap_G
     &                      / (25.D0 * gamma * t**2.D0)

c                 now, compute the density using its approximate form for
c                 small r/cap_R
                  cap_G = fudge_rho*rho2/rho1*
     &                      (r/cap_R)**(3.D0/(gamma-1.D0))

              end if

c             given dimensionless V and rho, compute their
c             physical values:

              velocity = 2.D0*r*cap_V/(5.D0*t)
              density = rho1 * cap_G

              end if

              rho_out (n) = density
              v_out   (n) = velocity
              p_out   (n) = pressure
        end do 

      end

      FUNCTION  func_TS (V, xi, gamma, nu1, nu2)

*//////////////////////////////////////////////////////////////////////
*/
*.  For the Taylor Sedov explosion problem.
*/  func is obtained from eq. (106.10) of Landau and Lifshitz,
*/  Fluid Mechanics, 2nd. edition.  Its root gives us the
*/  dimensionless velocity for a given position
*/
*//////////////////////////////////////////////////////////////////////

      implicit none

      real*8 func_TS
      real*8 V, xi, gamma, nu1, nu2

c      nu1 = -1.0*(13*gamma**2 - 7*gamma + 12)
c     &        /((3*gamma-1)*(2*gamma+1))
c      nu2 = 5*(gamma-1)/(2*gamma+1)

      func_TS = 1/(0.5*(gamma+1)*V)**2
     &         * ((gamma+1)/(7-gamma)*(5-(3*gamma-1)*V))**nu1
     &         * ((gamma+1)/(gamma-1)*(gamma*V-1))**nu2
     &         - xi**5.0

      return

      end


      FUNCTION rtbis_TS (V_low, V_high, xi, gamma, nu1, nu2, xacc)

*//////////////////////////////////////////////////////////////////////
*/
*/ This solves for the root of L&L (106.10).
*/ rtbis is borrowed from Numerical Recipes. It is a bisection 
*/ algorithm, which we use to solve for pm using a call to func.
*/
*/ Note that the arguments to rtbis have been altered and the value of
*/ JMAX increased. Otherwise, it is identical to the NR version.
*/
*//////////////////////////////////////////////////////////////////////

      implicit none

      INTEGER JMAX
      REAL*8 rtbis_TS,x1,x2,xacc,func_TS
      real*8 V_low, V_high, xi, gamma, nu1, nu2
      EXTERNAL func_TS
      PARAMETER (JMAX=300)
      INTEGER j
      REAL*8 dx,f,fmid,xmid

      x1 = V_low
      x2 = V_high

      fmid=func_TS(x2, xi, gamma, nu1, nu2)
      f=func_TS(x1, xi, gamma, nu1, nu2)
      if(f*fmid.ge.0.) pause 'root must be bracketed in rtbis'
      if(f.lt.0.)then
        rtbis_TS=x1
        dx=x2-x1
      else
        rtbis_TS=x2
        dx=x1-x2
      endif
      do j=1,JMAX
        dx=dx*5.D-1
        xmid=rtbis_TS+dx
        fmid=func_TS(xmid, xi, gamma, nu1, nu2)
        if(fmid.le.0.)rtbis_TS=xmid
        if(dabs(dx).lt.xacc .or. fmid.eq.0.) return
      enddo
      pause 'too many bisections in rtbis'
      END
