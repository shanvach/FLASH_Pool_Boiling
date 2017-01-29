! Potential flow around a cylinder
! Shizhao Wang
! March 8, 2015
! ===================================================

        subroutine sm_setPotFlowEllp(x,y,u,v,p)
        implicit none
        real :: x, y, u, v, p
        real :: u0, alpha, a1, b1, a, c
        real ::  xTmp, yTmp
        complex :: z, uz, dwdz, alphaz, zTmp

        u0 = 1.0
        alpha = acos(-1.0)/4.0
        !a1 = 1.005/2.0
        a1 = 0.5
        !b1 = 0.1007/2.0
        b1 = 0.25
        a  = 0.5*(a1+b1)
        c = sqrt(a1*a1-b1*b1)/2.0

        if(sqrt(x*x+y*y) < 1.0d-12) then
          x = 1.0d-12
          y = 1.0d-12
        endif

        alphaz = cmplx(cos(alpha), sin(alpha))

        z = cmplx(x,y)*alphaz
        zTmp = sqrt(z*z-4.0*c*c)
        if(real(z) < 0.0) then
          xTmp = real(zTmp)
          yTmp = imag(zTmp)
          zTmp = cmplx(-xTmp, -yTmp)
        endif
        dwdz = (u0*conjg(alphaz) - 4.0*a*a*alphaz/((z+zTmp)*(z+zTmp)))*(0.5+0.5*z/(zTmp))*alphaz
        u = real(dwdz)
        v = -imag(dwdz)
        p = -(u*u+v*v)

        return
        endSubroutine sm_setPotFlowEllp
