

! NOTE -- this was moved to ShockCyl from flash/source/utils/tools because the
!   only simulation using it was ShockCyl.

MODULE sim_ranluxModule

!     Subtract-and-borrow random number generator proposed by
!     Marsaglia and Zaman, implemented by F. James with the name
!     RCARRY in 1991, and later improved by Martin Luescher
!     in 1993 to produce "Luxury Pseudorandom Numbers".
!     Fortran 77 coded by F. James, 1993

!  References:
!  M. Luscher, Computer Physics Communications  79 (1994) 100
!  F. James, Computer Physics Communications 79 (1994) 111

!   LUXURY LEVELS.
!   ------ ------      The available luxury levels are:

!  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
!           and Zaman, very long period, but fails many tests.
!  level 1  (p=48): considerable improvement in quality over level 0,
!           now passes the gap test, but still fails spectral test.
!  level 2  (p=97): passes all known tests, but theoretically still
!           defective.
!  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
!           correlations have very small chance of being observed.
!  level 4  (p=389): highest possible luxury, all 24 bits chaotic.

!!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!!  Calling sequences for RANLUX:                                  ++
!!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
!!!!                   32-bit random floating point numbers between  ++
!!!!                   zero (not included) and one (also not incl.). ++
!!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
!!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
!!!!               which is integer between zero and MAXLEV, or if   ++
!!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
!!!!               should be set to zero unless restarting at a break++
!!!!               point given by output of RLUXAT (see RLUXAT).     ++
!!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
!!!!               which can be used to restart the RANLUX generator ++
!!!!               at the current point by calling RLUXGO.  K1 and K2++
!!!!               specify how many numbers were generated since the ++
!!!!               initialization with LUX and INT.  The restarting  ++
!!!!               skips over  K1+K2*E9   numbers, so it can be long.++
!!!!   A more efficient but less convenient way of restarting is by: ++
!!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
!!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
!!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
!!!!                 32-bit integer seeds, to be used for restarting ++
!!!!      ISVEC must be dimensioned 25 in the calling program        ++
!!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

IMPLICIT NONE

INTEGER            :: iseeds(24), isdext(25)
INTEGER, PARAMETER :: maxlev = 4, lxdflt = 3
INTEGER            :: ndskip(0:maxlev) = (/ 0, 24, 73, 199, 365 /)
INTEGER            :: next(24), igiga = 1000000000, jsdflt = 314159265
REAL, PARAMETER    :: twop12 = 4096.
INTEGER, PARAMETER :: itwo24 = 2**24, icons = 2147483563
INTEGER            :: luxlev = lxdflt, nskip, inseed, jseed
LOGICAL            :: notyet = .true.
INTEGER            :: in24 = 0, kount = 0, mkount = 0, i24 = 24, j24 = 10
REAL               :: seeds(24), carry = 0., twom24, twom12

!                            default
!  Luxury Level     0   1   2  *3*    4
!    ndskip        /0, 24, 73, 199, 365/
! Corresponds to p=24  48  97  223  389
!     time factor   1   2   3    6   10   on slow workstation
!                   1 1.5   2    3    5   on fast mainframe

PUBLIC notyet, i24, j24, carry, seeds, twom24, twom12, luxlev
PUBLIC nskip, ndskip, in24, next, kount, mkount, inseed


CONTAINS

!! --------------------------------------------------------------------------------


  SUBROUTINE sim_ranlux(rvec, lenv)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: lenv
    REAL, INTENT(OUT)   :: rvec(lenv)

    !     Local variables

    INTEGER             :: i, k, lp, ivec, isk
    REAL                :: uni

    !  NOTYET is .TRUE. if no initialization has been performed yet.
    !              Default Initialization by Multiplicative Congruential

    IF (notyet) THEN
       notyet = .false.
       jseed = jsdflt
       inseed = jseed
       WRITE (6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ', jseed
       luxlev = lxdflt
       nskip = ndskip(luxlev)
       lp = nskip + 24
       in24 = 0
       kount = 0
       mkount = 0
       WRITE (6,'(A,I2,A,I4)') ' RANLUX DEFAULT LUXURY LEVEL =  ', luxlev,   &
            '    p =', lp
       twom24 = 1.e0
       DO i = 1, 24
          twom24 = twom24 * 0.5e0
          k = jseed / 53668
          jseed = 40014 * (jseed-k*53668) - k * 12211
          IF (jseed < 0) jseed = jseed + icons
          iseeds(i) = MOD(jseed,itwo24)
       END DO
       twom12 = twom24 * 4096.e0
       DO i = 1, 24
          seeds(i) = REAL(iseeds(i)) * twom24
          next(i) = i - 1
       END DO
       next(1) = 24
       i24 = 24
       j24 = 10
       carry = 0.e0
       IF (seeds(24) == 0.e0) carry = twom24
    END IF

    !          The Generator proper: "Subtract-with-borrow",
    !          as proposed by Marsaglia and Zaman,
    !          Florida State University, March, 1989

    DO ivec = 1, lenv
       uni = seeds(j24) - seeds(i24) - carry
       IF (uni < 0.e0) THEN
          uni = uni + 1.e0
          carry = twom24
       ELSE
          carry = 0.e0
       END IF
       seeds(i24) = uni
       i24 = next(i24)
       j24 = next(j24)
       rvec(ivec) = uni
       !  small numbers (with less than 12 "significant" bits) are "padded".
       IF (uni < twom12) THEN
          rvec(ivec) = rvec(ivec) + twom24 * seeds(j24)
          !        and zero is forbidden in case someone takes a logarithm
          IF (rvec(ivec) == 0.e0) rvec(ivec) = twom24 * twom24
       END IF
       !        Skipping to luxury.  As proposed by Martin Luscher.
       in24 = in24 + 1
       IF (in24 == 24) THEN
          in24 = 0
          kount = kount + nskip
          DO isk = 1, nskip
             uni = seeds(j24) - seeds(i24) - carry
             IF (uni < 0.e0) THEN
                uni = uni + 1.0
                carry = twom24
             ELSE
                carry = 0.e0
             END IF
             seeds(i24) = uni
             i24 = next(i24)
             j24 = next(j24)
          END DO
       END IF
    END DO
    kount = kount + lenv
    IF (kount >= igiga) THEN
       mkount = mkount + 1
       kount = kount - igiga
    END IF
    RETURN

  END SUBROUTINE sim_ranlux


  SUBROUTINE sim_rluxin
    !           Subroutine to input and float integer seeds from previous run
    !     the following IF BLOCK added by Phillip Helbig, based on conversation
    !     with Fred James; an equivalent correction has been published by James.

    IMPLICIT NONE

    !     Local variables

    INTEGER             :: i, isd

    IF (notyet) THEN
       WRITE (6,'(A)') ' Proper results ONLY with initialisation from 25 ',  &
            'integers obtained with RLUXUT'
       notyet = .false.
    END IF

    twom24 = 1.e0
    DO i = 1, 24
       next(i) = i - 1
       twom24 = twom24 * 0.5e0
    END DO
    next(1) = 24
    twom12 = twom24 * 4096.e0
    WRITE (6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
    WRITE (6,'(5X,5I12)') isdext
    DO i = 1, 24
       seeds(i) = REAL(isdext(i)) * twom24
    END DO
    carry = 0.e0
    IF (isdext(25) < 0) carry = twom24
    isd = IABS(isdext(25))
    i24 = MOD(isd,100)
    isd = isd / 100
    j24 = MOD(isd,100)
    isd = isd / 100
    in24 = MOD(isd,100)
    isd = isd / 100
    luxlev = isd
    IF (luxlev <= maxlev) THEN
       nskip = ndskip(luxlev)
       WRITE (6,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ', luxlev
    ELSE IF (luxlev >= 24) THEN
       nskip = luxlev - 24
       WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:', luxlev
    ELSE
       nskip = ndskip(maxlev)
       WRITE (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ', luxlev
       luxlev = maxlev
    END IF
    inseed = -1
    RETURN

  END SUBROUTINE sim_rluxin

!! --------------------------------------------------------------------------------

  SUBROUTINE sim_rluxut
    !                    Subroutine to ouput seeds as integers

    IMPLICIT NONE

    !     Local variables

    INTEGER             :: i

    DO i = 1, 24
       isdext(i) = INT(seeds(i)*twop12*twop12)
    END DO
    isdext(25) = i24 + 100 * j24 + 10000 * in24 + 1000000 * luxlev
    IF (carry > 0.e0) isdext(25) = -isdext(25)
    RETURN

  END SUBROUTINE sim_rluxut

!! --------------------------------------------------------------------------------

  SUBROUTINE sim_rluxat(lout, inout, k1, k2)
    !                    Subroutine to output the "convenient" restart point

    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: lout, inout, k1, k2

    lout = luxlev
    inout = inseed
    k1 = kount
    k2 = mkount
    RETURN

  END SUBROUTINE sim_rluxat

!! --------------------------------------------------------------------------------

  SUBROUTINE sim_rluxgo(lux, ins, k1, k2, lwrite)
    !                    Subroutine to initialize from one or three integers

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: lux, ins, k1, k2
    LOGICAL, INTENT(IN) :: lwrite

    !     Local variables

    INTEGER             :: ilx, i, iouter, isk, k, inner, izip, izip2
    REAL                :: uni

    IF (lux < 0) THEN
       luxlev = lxdflt
    ELSE IF (lux <= maxlev) THEN
       luxlev = lux
    ELSE IF (lux < 24 .OR. lux > 2000) THEN
       luxlev = maxlev
       if ( lwrite ) &
            WRITE (6,'(A,I7)') ' [RLUXGO] Warning: Illegal luxury level : ', lux
    ELSE
       luxlev = lux
       DO ilx = 0, maxlev
          IF (lux == ndskip(ilx)+24) luxlev = ilx
       END DO
    END IF
    IF (luxlev <= maxlev) THEN
       nskip = ndskip(luxlev)
       if ( lwrite ) &
            WRITE (6,'(A,I2,A,I4)') ' [RLUXGO] Luxury level : ',luxlev,',  p=', nskip + 24
    ELSE
       nskip = luxlev - 24
       if ( lwrite ) &
            WRITE (6,'(A,I5)') ' [RLUXGO] Luxury level : ', luxlev
    END IF
    in24 = 0
    IF (ins < 0 .and. lwrite ) &
         WRITE (6,'(A)') &
         ' [RLUXGO] Warning: Illegal initialization, negative input seed'
    IF (ins >  0) THEN
       jseed = ins
       if ( lwrite ) &
            WRITE (6,'(A,3I12)') ' [RLUXG0] Initialization from seeds ', jseed, k1, k2
    ELSE
       jseed = jsdflt
       if ( lwrite ) &
            WRITE (6,'(A)') ' [RLUXGO] Initialized from default seed ',jseed
    END IF
    inseed = jseed
    notyet = .false.
    twom24 = 1.e0
    DO i = 1, 24
       twom24 = twom24 * 0.5e0
       k = jseed / 53668
       jseed = 40014 * (jseed-k*53668) - k * 12211
       IF (jseed < 0) jseed = jseed + icons
       iseeds(i) = MOD(jseed,itwo24)
    END DO
    twom12 = twom24 * 4096.e0
    DO i = 1, 24
       seeds(i) = REAL(iseeds(i)) * twom24
       next(i) = i - 1
    END DO
    next(1) = 24
    i24 = 24
    j24 = 10
    carry = 0.e0
    IF (seeds(24) == 0.e0) carry = twom24
    !        If restarting at a break point, skip K1 + IGIGA*K2
    !        Note that this is the number of numbers delivered to
    !        the user PLUS the number skipped (if luxury >  0).
    kount = k1
    mkount = k2
    IF (k1+k2 /= 0) THEN
       DO iouter = 1, k2 + 1
          inner = igiga
          IF (iouter == k2+1) inner = k1
          DO isk = 1, inner
             uni = seeds(j24) - seeds(i24) - carry
             IF (uni < 0.e0) THEN
                uni = uni + 1.e0
                carry = twom24
             ELSE
                carry = 0.e0
             END IF
             seeds(i24) = uni
             i24 = next(i24)
             j24 = next(j24)
          END DO
       END DO
       !         Get the right value of IN24 by direct calculation
       in24 = MOD(kount,nskip+24)
       IF (mkount > 0) THEN
          izip = MOD(igiga, nskip+24)
          izip2 = mkount * izip + in24
          in24 = MOD(izip2, nskip+24)
       END IF
       !       Now IN24 had better be between zero and 23 inclusive
       IF (in24 > 23) THEN
          if ( lwrite ) &
               WRITE (6,'(A/A,3I11,A,I5)') &
               ' [RLUXGO] ERROR: while restarting with :', '  The values', ins, &
               k1, k2, ' cannot occur at luxury level', luxlev
          in24 = 0
       END IF
    END IF
    RETURN

  END SUBROUTINE sim_rluxgo

END MODULE sim_ranluxModule
