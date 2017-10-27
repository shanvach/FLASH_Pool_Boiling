subroutine Heat_getQmicro(qmic,fmic,dxmin)

     use Heat_AD_data, only: ht_Pr, ht_St, ht_Ab, ht_Bb, ht_Cb, ht_Twall_low, ht_Tsat

     use IncompNS_data, only: ins_invRe

     use Multiphase_data, only: mph_rho1, mph_rho2, mph_sten 

     implicit none

     real, intent(inout) :: qmic,fmic
     real, intent(in)    :: dxmin

     real :: Re, Pr, St, rho, We, Pe, Ab, Bb, Cb, Tw, Ts
     real :: dr,r,psi,pi,step
     integer :: N,i
     real, allocatable, dimension(:) :: z1,z2,z3,z4,q

     qmic = 0.0
     fmic = 0.0

     pi  = acos(-1.d0)

     psi = 35*pi/180

     Re  = 1.0/ins_invRe
     Pr  = ht_Pr
     We  = 1.0/mph_sten
     rho = mph_rho1/mph_rho2
     St  = ht_St
     Ab  = ht_Ab
     Bb  = ht_Bb
     Cb  = ht_Cb
     Tw  = ht_Twall_low
     Ts  = ht_Tsat

     r  = dxmin/(2.0*tan(psi))
     !dr = 0.95d-4
     !dr = 0.87d-4
     !dr = 0.40d-4
     dr = 0.7d-4

     step = dr
     N = r/dr

     allocate(z1(N))
     allocate(z2(N))
     allocate(z3(N))
     allocate(z4(N))
     allocate(q(N))

     z1(1) = dxmin/2.0
     z2(1) = tan(psi)
     z3(1) = Ab/((dxmin/2.0)**3)
     z4(1) = (Re*Pr*z2(1)*Ab)/(St*z1(1))
     dr    = -dr
     q(1)  = (Tw-Ts-(Bb/Re)*z3(1))/(z1(1) + Cb/rho)
     qmic = qmic - dr*q(1)

     do i=2,N

        z1(i) = z1(i-1) + dr*z2(i-1)
        z2(i) = z2(i-1) + dr*((z3(i-1) - Ab/(z1(i-1)**3))/(Re/We))*((1+z2(i-1)**2)**(3.0/2.0))
        z3(i) = z3(i-1) + dr*(3.0*St/(Re*Pr))*(z4(i-1)/(z1(i-1)**3))
        z4(i) = z4(i-1) + dr*(Ts - Tw + (Bb/Re)*z3(i-1))/(z1(i-1) + Cb/rho)

        q(i)  = (Tw-Ts-(Bb/Re)*z3(i))/(z1(i) + Cb/rho)
        qmic = qmic - step*q(i)
        fmic = fmic + step*(z3(i)-(Ab/(z1(i)**3)))/(Re/We)

     end do

     deallocate(z1,z2,z3,z4,q)

end subroutine Heat_getQmicro
