subroutine bsetup(myPE)

  use Simulation_data, ONLY: sim_Bx, sim_By, sim_Bz, sim_Bdx, &
       sim_Bdy, sim_Bdz, sim_nBzones, lbox, niq, sim_Bxcoord, &
       sim_Bycoord, sim_Bzcoord, sim_ptdirn, &
       sim_Bmag, sim_pressureNormalize, sim_plasmaBeta, &
       sim_xMin, sim_yMin, sim_zMin, r1, pres1, numPoints1
  
  IMPLICIT NONE

#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: myPE

  integer, parameter :: NN = 256

  real, dimension(NN,NN,NN) :: Bx
  real, dimension(NN,NN,NN) :: By
  real, dimension(NN,NN,NN) :: Bz
  complex, dimension(NN/2,NN,NN) :: Bxc
  complex, dimension(NN/2,NN,NN) :: Byc
  complex, dimension(NN/2,NN,NN) :: Bzc
  COMPLEX, dimension(NN,NN) :: Bxspeq
  COMPLEX, dimension(NN,NN) :: Byspeq
  COMPLEX, dimension(NN,NN) :: Bzspeq
  
  complex :: a, b, c, kr(3), Bk
  INTEGER :: i,j,k,n,ii,jj,kk
  INTEGER :: idum    
  COMPLEX, external :: getBc
  REAL :: modk, Bavg,Bmag,xx,yy,zz,pres,rr,dr,totpres
  real :: totmagp, beta
  REAL, external :: getmodk
  real, external :: interpolate

  equivalence(Bx,Bxc)
  equivalence(By,Byc)
  equivalence(Bz,Bzc)

  Bx = 0.0
  By = 0.0
  Bz = 0.0

  niq = NN/2+1

  idum = -1000000+time()

  print *, "Begin setting up B-field"

  dr = sqrt(sim_Bdx*sim_Bdx+sim_Bdy*sim_Bdy+sim_Bdz*sim_Bdz)

!------------------------------------------------------------------------------
! Setup of the FT (Bc) of the B-field
!
!
! Because B is real Bc is Hermitian: Bc(k1,k2,k3) = Bc*(-k1,-k2,-k3).
! For kN=k_Niquist it holds Bc(kN,k2,k3) = Bc*(kN,-k2,-k3),
! as well as Bc(0,k2,k3) = Bc*(0,-k2,-k3).
! 
!  k ---------------------------------------    i=1 or i =niq
!    |   |              |   |              |
!    |   |              |   |              |
!    | 2 |       5      | 2 |       4      |
!    |   |              |   |              |
!    |   |              |   |              |
!    |   |              |   |              |
!    ---------------------------------------
!    | 1 |       3      | 1 |       3      |
!    ---------------------------------------
!    |   |              |   |              |
!    |   |              |   |              |
!    | 2 |       4      | 2 |       5      |
!    |   |              |   |              |
!    |   |              |   |              |
!    |   |              |   |              |
!    ---------------------------------------
!    | 1 |       3      | 1 |       3      |
!    ---------------------------------------
!                                          j

! --- region 1 above
  i=1
  DO j=1,niq,niq-1
     DO k=1,niq,niq-1
        modk = getmodk(i,j,k)
        Bx(i,j,k) = getBc(modk,1,idum)
        By(i,j,k) = getBc(modk,1,idum)
        Bz(i,j,k) = getBc(modk,1,idum)
     ENDDO
  ENDDO

  i=niq
  DO j=1,niq,niq-1
     DO k=1,niq,niq-1
        modk = getmodk(i,j,k)
        Bxspeq(j,k)=getBc(modk,1,idum)
        Byspeq(j,k)=getBc(modk,1,idum)
        Bzspeq(j,k)=getBc(modk,1,idum)
     ENDDO
  ENDDO
  
! --- region 2 above
  i=1
  DO j=1,niq,niq-1
     DO k=2,niq-1
        modk = getmodk(i,j,k)
        Bxc(i,j,k)=getBc(modk,0,idum)
        Bxc(i,j,NN-k+2)=CONJG(Bxc(i,j,k))
        Byc(i,j,k)=getBc(modk,0,idum)
        Byc(i,j,NN-k+2)=CONJG(Byc(i,j,k))
        Bzc(i,j,k)=getBc(modk,0,idum)
        Bzc(i,j,NN-k+2)=CONJG(Bzc(i,j,k))
     ENDDO
  ENDDO
  
  i=niq
  DO j=1,niq,niq-1
     DO k=2,niq-1
        modk = getmodk(i,j,k)
        Bxspeq(j,k)=getBc(modk,0,idum)
        Bxspeq(j,NN-k+2)=CONJG(Bxspeq(j,k))
        Byspeq(j,k)=getBc(modk,0,idum)
        Byspeq(j,NN-k+2)=CONJG(Byspeq(j,k))
        Bzspeq(j,k)=getBc(modk,0,idum)
        Bzspeq(j,NN-k+2)=CONJG(Bzspeq(j,k))
     ENDDO
 ENDDO
 
! --- region 3 above
 i=1
 DO j=2,niq-1
    DO k=1,niq,niq-1
       modk = getmodk(i,j,k)
       Bxc(i,j,k)=getBc(modk,0,idum)
       Bxc(i,NN-j+2,k)=CONJG(Bxc(i,j,k))
       Byc(i,j,k)=getBc(modk,0,idum)
       Byc(i,NN-j+2,k)=CONJG(Byc(i,j,k))
       Bzc(i,j,k)=getBc(modk,0,idum)
       Bzc(i,NN-j+2,k)=CONJG(Bzc(i,j,k))
    ENDDO
 ENDDO
 
 i=niq
 DO j=2,niq-1
    DO k=1,niq,niq-1
       modk = getmodk(i,j,k)
       Bxspeq(j,k)=getBc(modk,0,idum)
       Bxspeq(NN-j+2,k)=CONJG(Bxspeq(j,k))
       Byspeq(j,k)=getBc(modk,0,idum)
       Byspeq(NN-j+2,k)=CONJG(Byspeq(j,k))
       Bzspeq(j,k)=getBc(modk,0,idum)
       Bzspeq(NN-j+2,k)=CONJG(Bzspeq(j,k))
    ENDDO
 ENDDO
    
! --- region 4 above
 i=1
 DO j=2,niq-1
    DO k=2,niq-1
       modk = getmodk(i,j,k)
       Bxc(i,j,k)=getBc(modk,0,idum)
       Bxc(i,NN-j+2,NN-k+2)=CONJG(Bxc(i,j,k))
       Byc(i,j,k)=getBc(modk,0,idum)
       Byc(i,NN-j+2,NN-k+2)=CONJG(Byc(i,j,k))
       Bzc(i,j,k)=getBc(modk,0,idum)
       Bzc(i,NN-j+2,NN-k+2)=CONJG(Bzc(i,j,k))
    ENDDO
 ENDDO
      
 i=niq
 DO j=2,niq-1
    DO k=2,niq-1
       modk = getmodk(i,j,k)
       Bxspeq(j,k)=getBc(modk,0,idum)
       Bxspeq(NN-j+2,NN-k+2)=CONJG(Bxspeq(j,k))
       Byspeq(j,k)=getBc(modk,0,idum)
       Byspeq(NN-j+2,NN-k+2)=CONJG(Byspeq(j,k))
       Bzspeq(j,k)=getBc(modk,0,idum)
       Bzspeq(NN-j+2,NN-k+2)=CONJG(Bzspeq(j,k))
    ENDDO
 ENDDO

! --- region 5 above
 i=1
 DO j=niq+1,NN
    DO k=2,niq-1
       modk = getmodk(i,j,k)
       Bxc(i,j,k)=getBc(modk,0,idum)
       Bxc(i,NN-j+2,NN-k+2)=CONJG(Bxc(i,j,k))
       Byc(i,j,k)=getBc(modk,0,idum)
       Byc(i,NN-j+2,NN-k+2)=CONJG(Byc(i,j,k))
       Bzc(i,j,k)=getBc(modk,0,idum)
       Bzc(i,NN-j+2,NN-k+2)=CONJG(Bzc(i,j,k))
    ENDDO
 ENDDO

 i=niq
 DO j=niq+1,NN
    DO k=2,niq-1
       modk = getmodk(i,j,k)
       Bxspeq(j,k)=getBc(modk,0,idum)
       Bxspeq(NN-j+2,NN-k+2)=CONJG(Bxspeq(j,k))
       Byspeq(j,k)=getBc(modk,0,idum)
       Byspeq(NN-j+2,NN-k+2)=CONJG(Byspeq(j,k))
       Bzspeq(j,k)=getBc(modk,0,idum)
       Bzspeq(NN-j+2,NN-k+2)=CONJG(Bzspeq(j,k))
    ENDDO
 ENDDO

! --- cells without Hermitian counterpart
 DO i=2,niq-1
    DO j=1,NN
       DO k=1,NN
          modk = getmodk(i,j,k)
          Bxc(i,j,k)=getBc(modk,0,idum)
          Byc(i,j,k)=getBc(modk,0,idum)
          Bzc(i,j,k)=getBc(modk,0,idum)
       ENDDO
    ENDDO
 ENDDO

!------------------------------------------------------------------------------
! FT of Bc -> B

 CALL rlft3(Bx,Bxspeq,NN,NN,NN,-1)
 CALL rlft3(By,Byspeq,NN,NN,NN,-1)
 CALL rlft3(Bz,Bzspeq,NN,NN,NN,-1)

!------------------------------------------------------------------------------
! Normalize 

 Bavg = 0.0
 do i = 1, sim_nBzones
    do j = 1, sim_nBzones
       do k = 1, sim_nBzones
          
          Bavg = Bavg + sqrt(Bx(i,j,k)**2+By(i,j,k)**2+Bz(i,j,k)**2)
          
       enddo
    enddo
 enddo
 
 Bavg = Bavg/(sim_nBzones**3)
 
 Bx(:,:,:) = (Bx(:,:,:) / Bavg)*sim_Bmag
 By(:,:,:) = (By(:,:,:) / Bavg)*sim_Bmag
 Bz(:,:,:) = (Bz(:,:,:) / Bavg)*sim_Bmag
 
!------------------------------------------------------------------------------
! FT of B -> Bc

 CALL rlft3(Bx,Bxspeq,NN,NN,NN,1)
 CALL rlft3(By,Byspeq,NN,NN,NN,1)
 CALL rlft3(Bz,Bzspeq,NN,NN,NN,1)

!------------------------------------------------------------------------------

! Divergence cleaning

 do i = 1, NN/2
    do j = 1, NN
       do k = 1, NN

          a = Bxc(i,j,k)
          b = Byc(i,j,k)
          c = Bzc(i,j,k)

          call getks(i,j,k,kr)

          Bk = conjg(kr(1))*a + conjg(kr(2))*b + conjg(kr(3))*c

          Bxc(i,j,k) = Bxc(i,j,k) - kr(1)*Bk
          Byc(i,j,k) = Byc(i,j,k) - kr(2)*Bk
          Bzc(i,j,k) = Bzc(i,j,k) - kr(3)*Bk

       enddo
    enddo
 enddo

 i = niq

 do j = 1, NN
    do k = 1, NN

       a = Bxspeq(j,k)
       b = Byspeq(j,k)
       c = Bzspeq(j,k)

       call getks(i,j,k,kr)

       Bk = conjg(kr(1))*a + conjg(kr(2))*b + conjg(kr(3))*c

       Bxspeq(j,k) = Bxspeq(j,k) - kr(1)*Bk
       Byspeq(j,k) = Byspeq(j,k) - kr(2)*Bk
       Bzspeq(j,k) = Bzspeq(j,k) - kr(3)*Bk

    enddo
 enddo

!------------------------------------------------------------------------------
! FT of Bc -> B

 CALL rlft3(Bx,Bxspeq,NN,NN,NN,-1)
 CALL rlft3(By,Byspeq,NN,NN,NN,-1)
 CALL rlft3(Bz,Bzspeq,NN,NN,NN,-1)

 Bx = Bx * (2./NN**3)
 By = By * (2./NN**3)
 Bz = Bz * (2./NN**3)

 sim_Bx = 0.0
 sim_By = 0.0
 sim_Bz = 0.0

 if (sim_pressureNormalize) then

    totmagp = sum(Bx*Bx+By*By+Bz*Bz)/2.
    
    beta = totpres/totmagp
    
    Bx = Bx / sqrt(sim_plasmaBeta/beta) 
    By = By / sqrt(sim_plasmaBeta/beta)
    Bz = Bz / sqrt(sim_plasmaBeta/beta)
    
 endif

 if (NN >= sim_nBzones) then

    sim_Bx = Bx(1:sim_nBzones,1:sim_nBzones,1:sim_nBzones)
    sim_By = By(1:sim_nBzones,1:sim_nBzones,1:sim_nBzones)
    sim_Bz = Bz(1:sim_nBzones,1:sim_nBzones,1:sim_nBzones)
    
 else

    sim_Bx(1:NN,1:NN,1:NN) = Bx
    sim_By(1:NN,1:NN,1:NN) = By
    sim_Bz(1:NN,1:NN,1:NN) = Bz

    sim_Bx(NN+1:sim_nBzones,1:NN,1:NN) = Bx(1:NGUARD,1:NN,1:NN) 
    sim_By(NN+1:sim_nBzones,1:NN,1:NN) = By(1:NGUARD,1:NN,1:NN)
    sim_Bz(NN+1:sim_nBzones,1:NN,1:NN) = Bz(1:NGUARD,1:NN,1:NN)

    sim_Bx(1:NN,NN+1:sim_nBzones,1:NN) = Bx(1:NN,1:NGUARD,1:NN)
    sim_By(1:NN,NN+1:sim_nBzones,1:NN) = By(1:NN,1:NGUARD,1:NN)
    sim_Bz(1:NN,NN+1:sim_nBzones,1:NN) = Bz(1:NN,1:NGUARD,1:NN)

    sim_Bx(1:NN,1:NN,NN+1:sim_nBzones) = Bx(1:NN,1:NN,1:NGUARD)
    sim_By(1:NN,1:NN,NN+1:sim_nBzones) = By(1:NN,1:NN,1:NGUARD)
    sim_Bz(1:NN,1:NN,NN+1:sim_nBzones) = Bz(1:NN,1:NN,1:NGUARD)

    sim_Bx(NN+1:sim_nBzones,NN+1:sim_nBzones,1:NN) = Bx(1:NGUARD,1:NGUARD,1:NN)
    sim_By(NN+1:sim_nBzones,NN+1:sim_nBzones,1:NN) = By(1:NGUARD,1:NGUARD,1:NN)
    sim_Bz(NN+1:sim_nBzones,NN+1:sim_nBzones,1:NN) = Bz(1:NGUARD,1:NGUARD,1:NN)

    sim_Bx(1:NN,NN+1:sim_nBzones,NN+1:sim_nBzones) = Bx(1:NN,1:NGUARD,1:NGUARD)
    sim_By(1:NN,NN+1:sim_nBzones,NN+1:sim_nBzones) = By(1:NN,1:NGUARD,1:NGUARD)
    sim_Bz(1:NN,NN+1:sim_nBzones,NN+1:sim_nBzones) = Bz(1:NN,1:NGUARD,1:NGUARD)

    sim_Bx(NN+1:sim_nBzones,1:NN,NN+1:sim_nBzones) = Bx(1:NGUARD,1:NN,1:NGUARD)
    sim_By(NN+1:sim_nBzones,1:NN,NN+1:sim_nBzones) = By(1:NGUARD,1:NN,1:NGUARD)
    sim_Bz(NN+1:sim_nBzones,1:NN,NN+1:sim_nBzones) = Bz(1:NGUARD,1:NN,1:NGUARD)

    sim_Bx(NN+1:sim_nBzones,NN+1:sim_nBzones,NN+1:sim_nBzones) = Bx(1:NGUARD,1:NGUARD,1:NGUARD)
    sim_By(NN+1:sim_nBzones,NN+1:sim_nBzones,NN+1:sim_nBzones) = By(1:NGUARD,1:NGUARD,1:NGUARD)
    sim_Bz(NN+1:sim_nBzones,NN+1:sim_nBzones,NN+1:sim_nBzones) = Bz(1:NGUARD,1:NGUARD,1:NGUARD)


 endif
 
!------------------------------------------------------------------------------
 
 return

END subroutine bsetup 

COMPLEX FUNCTION getBc(modk,rn,idum)

  use Simulation_data, ONLY: sim_lMin, sim_lMax

  implicit none
    
  integer, parameter :: NN = 256
  REAL, intent(in) ::    modk
  INTEGER, intent(in) :: rn       ! rn =1 implies Bc real => phi=0./pi
  INTEGER, intent(in) ::  idum    
  REAL    sigma2
  REAL    a,phi
  REAL    ran3
  REAL    pi
  PARAMETER (pi = 3.4159265355)
  real :: k0
  real :: k1
  !parameter (k0 = 2.*pi/0.043)
  !parameter (k1 = 2.*pi/0.5)

  k0 = 2.*pi/(sim_lMin*1.0e-3)
  k1 = 2.*pi/(sim_lMax*1.0e-3)

  if (modk == 0.0) then
     sigma2 = 0.0
  else
     sigma2 = exp(-k1/modk)*modk**(-11./6.)*exp(-(modk/k0)**2)
  endif

  sigma2 = sigma2*sigma2

  a      = SQRT((-2.)*sigma2*LOG(ran3(idum)))
  phi    = 2.*pi*ran3(idum)
  IF (rn.EQ.1) THEN
     IF (phi.LT.pi) THEN
        getBc  = CMPLX(a,0.)
     ELSE
        getBc  = CMPLX(a*(-1),0)
     ENDIF
  ELSE
     getBc  = CMPLX(a*COS(phi),a*SIN(phi))
  ENDIF
  
  RETURN
END FUNCTION getBc
      
!****************************************************************************
! getmodk returns the module of K
!
! i> niq => ki = -k(lbox+2-i)
REAL FUNCTION getmodk(i,j,k)

  use Simulation_data, ONLY: sim_Bdx, niq, lbox

  implicit none

#include "constants.h"
  
  INTEGER, intent(in) :: i,j,k
  INTEGER ii,jj,kk

  integer, parameter :: NN = 256

  real, save :: dx

  dx = sim_Bdx / 3.0856e24

  IF (i.LE.niq) THEN
     ii = i - 1
  ELSE
     ii = i - NN -1
  ENDIF

  IF (j.LE.niq) THEN
     jj = j - 1
  ELSE
     jj = j - NN -1
  ENDIF

  IF (k.LE.niq) THEN
     kk = k - 1
  ELSE
     kk = k - NN -1
  ENDIF

  getmodk = 2.*PI*SQRT(REAL(ii**2+jj**2+kk**2))/(NN*dx)
  
  RETURN
END FUNCTION getmodk

subroutine getks(i,j,k,kr)

  use Simulation_data, ONLY : lbox, niq

  implicit none

#include "constants.h"

  INTEGER, intent(IN) :: i,j,k
  INTEGER :: ii,jj,kk

  integer, parameter :: NN = 256

  real :: kmod, kreal, kimag

  complex, dimension(3), intent(OUT) :: kr

  IF (i.LE.niq) THEN
     ii = i - 1
  ELSE
     ii = i - NN -1
  ENDIF

  IF (j.LE.niq) THEN
     jj = j - 1
  ELSE
     jj = j - NN -1
  ENDIF

  IF (k.LE.niq) THEN
     kk = k - 1
  ELSE
     kk = k - NN -1
  ENDIF

  if (ii == 0 .and. jj == 0 .and. kk == 0) then
     kr(:) = cmplx(0.0,0.0)
  else
     kr(1) = cmplx(sin(2*PI*ii/real(NN)),1.-cos(2*PI*ii/real(NN)))
     kr(2) = cmplx(sin(2*PI*jj/real(NN)),1.-cos(2*PI*jj/real(NN)))
     kr(3) = cmplx(sin(2*PI*kk/real(NN)),1.-cos(2*PI*kk/real(NN)))
     kmod = sqrt(real(conjg(kr(1))*kr(1) + &
          conjg(kr(2))*kr(2) + &
          conjg(kr(3))*kr(3)))
     kr(:) = kr(:)/kmod
  endif

  return

end subroutine getks
