! Potential flow around a cylinder
! Shizhao Wang
! March 8, 2015
! ===================================================

        program flowEll
        implicit none
        real :: u0, alpha, a, b, c, u, v, p, x, y
        real :: x0, y0, dx, dy, e0, xTmp, yTmp
        integer :: nX, nY, i, j
        complex :: z, uz, dwdz, alphaz, zeta, zeta0, zTmp

        u0 = 1.0
        alpha = acos(-1.0)/4.0
        a = 0.5
        b = 0.25
        c = sqrt(a*a-b*b)/2.0

        alphaz = cmplx(cos(alpha), sin(alpha))

        x0 = -2.0+0.02
        y0 = -2.0+0.02
        dx = 0.04
        dy = 0.04
        nX = 100
        nY = 100

        open(100, file='flowEll.dat')
        write(100,*) 'VARIABLES="x", "y", "u", "v", "p"'
        write(100,*) 'ZONE, I=', nX, ', J=', nY, ', F=POINT'  
        do j = 1, nY
          y = y0 + (j-1)*dy
          do i = 1, nX
            x = x0 + (i-1)*dx
            z = cmplx(x,y)*alphaz
            xTmp = real(z*z-4.0*c*c)
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
            write(100, *) x, y, u, v, p
          enddo
        enddo
        close(100)

        contains
        ! ================================================
        function sh(z0)
        complex :: sh, z0
        real :: x, y
          x = real(z0)
          y = imag(z0)
          sh = cmplx(sinh(x)*cos(y),cosh(x)*sin(y))

        endfunction sh

        endProgram flowEll
