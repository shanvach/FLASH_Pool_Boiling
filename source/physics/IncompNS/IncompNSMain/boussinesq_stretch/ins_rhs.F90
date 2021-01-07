
      SUBROUTINE ins_rhs2d(uni,vni,tni,ru1,ix1,ix2,jy1,jy2,dx,dy,ru,rv)

     !***************************************************************
     ! This subroutine computes the discretization of the RHS of the 
     ! Helmholtz equation on a staggered uniform grid.
     !
     ! Input:  uni,vni     = velocity at timestep n
     !         tni         = temperature at timestep n
     !         ru1         = molecular viscosity
     !         ix1,ix2     = starting and ending x indices
     !         jy1,jy2     = starting and ending y indices
     !         dx,dy       = grid spacing in x and y directions
     !
     ! Output: ru,rv    = u and v momentum for Helmholtz RHS
     !**************************************************************

      use IncompNS_data, ONLY : ins_isCoupled

      implicit none

#include "constants.h"

      INTEGER, INTENT(IN) :: ix1, ix2, jy1, jy2
      REAL,    INTENT(IN) :: ru1
      REAL, DIMENSION(:,:),   INTENT(IN)  :: dx, dy
      REAL, DIMENSION(:,:,:), INTENT(IN)  :: uni, vni, tni
      REAL, DIMENSION(:,:,:), INTENT(OUT) :: ru, rv

      INTEGER:: i, j
      ! x-component variables
      REAL:: uxplus, uxminus, vxplus, vxminus, wxplus, wxminus, &
             uyplus, uyminus
      REAL:: dudxp, dudxm, dudyp, dudym, dvdxp, dvdxm
      REAL:: tvjp, tvjm
      REAL:: txxp, txxm, tyyp, tyym
      REAL:: txyp, txym
      
      ! new y-component variables
      REAL:: vyplus, vyminus
      REAL:: dvdyp, dvdym
      REAL:: tvip, tvim, ty, ty_m
      INTEGER, parameter :: kz1 = 1

      ! coupled term multiplier
      if (ins_isCoupled) then
        ty_m = 1.
      else
        ty_m = 0.
      endif

      !++++++++++  U-COMPONENT  ++++++++++
       do j = jy1,jy2
          do i = ix1,ix2+1
             ! get velocities at 1/2 locations
             uxplus = (uni(i+1,j,kz1) + uni(i,j,kz1))*0.5
             uxminus = (uni(i,j,kz1) + uni(i-1,j,kz1))*0.5

             vxplus = (vni(i,j+1,kz1) + vni(i-1,j+1,kz1))*0.5
             vxminus = (vni(i,j,kz1) + vni(i-1,j,kz1))*0.5

             uyplus = (uni(i,j+1,kz1) + uni(i,j,kz1))*0.5
             uyminus = (uni(i,j,kz1) + uni(i,j-1,kz1))*0.5

             ! get derivatives at 1/2 locations
             dudxp = (uni(i+1,j,kz1) - uni(i,j,kz1)) * dx(i,   CENTER)
             dudxm = (uni(i,j,kz1) - uni(i-1,j,kz1)) * dx(i-1, CENTER)
             dudyp = (uni(i,j+1,kz1) - uni(i,j,kz1)) * dy(j,   RIGHT_EDGE)
             dudym = (uni(i,j,kz1) - uni(i,j-1,kz1)) * dy(j,   LEFT_EDGE)
             !dvdxp = (vni(i,j+1,kz1) - vni(i-1,j+1,kz1)) * dx(i,   LEFT_EDGE)
             !dvdxm = (vni(i,j,kz1) - vni(i-1,j,kz1))     * dx(i,   LEFT_EDGE)

             ! flux of normal total stresses
             txxp = ru1*dudxp
             txxm = ru1*dudxm
             tyyp = ru1*dudyp
             tyym = ru1*dudym

             ! calculate RHS for u-momentum
             ru(i,j,kz1) =                                                    &
                          - (uxplus*uxplus - uxminus*uxminus)*dx(i,LEFT_EDGE) &! advection term
                          - (vxplus*uyplus - vxminus*uyminus)*dy(j,CENTER)    &                          
                          + (txxp - txxm)*dx(i,LEFT_EDGE)                     &! diffusion - normal terms 
                          + (tyyp - tyym)*dy(j,CENTER)
         enddo
       enddo

    !++++++++++  V-COMPONENT  ++++++++++

       do j = jy1,jy2+1
          do i = ix1,ix2
             ! get velocities at 1/2 locations
             vxplus = (vni(i+1,j,kz1) + vni(i,j,kz1))*0.5
             vxminus = (vni(i,j,kz1) + vni(i-1,j,kz1))*0.5

             vyplus = (vni(i,j+1,kz1) + vni(i,j,kz1))*0.5
             vyminus = (vni(i,j,kz1) + vni(i,j-1,kz1))*0.5

             uyplus = (uni(i+1,j,kz1) + uni(i+1,j-1,kz1))*0.5
             uyminus = (uni(i,j,kz1) + uni(i,j-1,kz1))*0.5

             ! get temperature at v locations
             ty = (tni(i,j,kz1) + tni(i,j-1,kz1))*0.5
             
             ! get derivatives at 1/2 locations
             dvdxp = (vni(i+1,j,kz1) - vni(i,j,kz1)) * dx(i,   RIGHT_EDGE)
             dvdxm = (vni(i,j,kz1) - vni(i-1,j,kz1)) * dx(i,   LEFT_EDGE)
             dvdyp = (vni(i,j+1,kz1) - vni(i,j,kz1)) * dy(j,   CENTER)
             dvdym = (vni(i,j,kz1) - vni(i,j-1,kz1)) * dy(j-1, CENTER)
             !dudyp = (uni(i+1,j,kz1) - uni(i+1,j-1,kz1)) * dy(j, LEFT_EDGE)
             !dudym = (uni(i,j,kz1) - uni(i,j-1,kz1))     * dy(j, LEFT_EDGE)  

             ! flux of normal total stresses
             txxp = ru1*dvdxp
             txxm = ru1*dvdxm
             tyyp = ru1*dvdyp
             tyym = ru1*dvdym

             ! calculate RHS for v-momentum
             rv(i,j,kz1) =                                                    &
                          - (uyplus*vxplus - uyminus*vxminus)*dx(i,CENTER)    &! advection term
                          - (vyplus*vyplus - vyminus*vyminus)*dy(j,LEFT_EDGE) &
                          + (txxp - txxm)*dx(i,CENTER)                        &! diffusion - normal terms
                          + (tyyp - tyym)*dy(j,LEFT_EDGE)                     &!
                          + ty * ty_m                                          ! energy momentum coupling
          enddo
       enddo

       END SUBROUTINE ins_rhs2d




  SUBROUTINE ins_rhs3d(uni,vni,wni,tni,tv,ru1,ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz,ru,rv,rw)

  !*****************************************************************
  ! This subroutine computes the centered discretization of the RHS 
  ! of the momentum equation (advection + viscous terms) on a 
  ! staggered uniform grid based on the Paramesh grid structure.
  !
  ! Input:  uni,vni,wni = velocity at timestep n
  !         tv          = eddy viscosity
  !         ru1         = molecular viscosity
  !         ix1,ix2     = starting and ending x indices
  !         jy1,jy2     = starting and ending y indices
  !         kz1,kz2     = starting and ending z indices
  !         dx,dy,dz    = grid spacing in x, y, and z directions
  !
  ! Output: ru,rv,rw    = RHS of u, v, and w momentum equations
  !
  ! E. Balaras   July 1999
  ! P. Rabenold  August 2006
  ! A. Lentner   December 2019
  !**************************************************************
    
  use IncompNS_data, ONLY : ins_isCoupled

  implicit none
  
  INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2, kz1, kz2
  REAL, INTENT(IN):: ru1
  REAL, DIMENSION(:,:),   INTENT(IN)  :: dx, dy, dz
  REAL, DIMENSION(:,:,:), INTENT(IN)  :: uni, vni, wni, tni, tv
  REAL, DIMENSION(:,:,:), INTENT(OUT) :: ru, rv, rw

  INTEGER:: i, j, k

  ! x-component variables
  REAL :: uxplus, uxminus, vxplus, vxminus, wxplus, wxminus, &
          uyplus, uyminus, uzplus, uzminus
  REAL :: dudxp, dudxm, dudyp, dudym, dudzp, dudzm, dvdxp, dvdxm, &
          dwdxp, dwdxm
  REAL :: tvjp, tvjm, tvkp, tvkm
  REAL :: txxp, txxm, tyyp, tyym, tzzp, tzzm
  REAL :: txyp, txym, txzp, txzm

  ! additional y-component variables
  REAL :: vyplus, vyminus, vzplus, vzminus, wyplus, wyminus
  REAL :: dvdyp, dvdym, dvdzp, dvdzm, dwdyp, dwdym
  REAL :: tvip, tvim
  REAL :: tyzp, tyzm

  ! additional z-component variables
  REAL :: wzplus, wzminus
  REAL :: dwdzp, dwdzm
  REAL :: tz, tz_m

  ! coupled term multiplier
  if (ins_isCoupled) then
    tz_m = 1.
  else
    tz_m = 0.
  endif

  !++++++++++  U-COMPONENT  ++++++++++
  do k = kz1,kz2
    do j = jy1,jy2
      do i = ix1,ix2+1

        ! get velocities at 1/2 locations
        uxplus  = (uni(i+1,j  ,k  ) + uni(i  ,j  ,k  ))*0.5
        uxminus = (uni(i  ,j  ,k  ) + uni(i-1,j  ,k  ))*0.5

        uyplus  = (uni(i  ,j+1,k  ) + uni(i  ,j  ,k  ))*0.5
        uyminus = (uni(i  ,j  ,k  ) + uni(i  ,j-1,k  ))*0.5

        uzplus  = (uni(i  ,j  ,k+1) + uni(i  ,j  ,k  ))*0.5
        uzminus = (uni(i  ,j  ,k  ) + uni(i  ,j  ,k-1))*0.5

        vxplus  = (vni(i  ,j+1,k  ) + vni(i-1,j+1,k  ))*0.5
        vxminus = (vni(i  ,j  ,k  ) + vni(i-1,j  ,k  ))*0.5

        wxplus  = (wni(i  ,j  ,k+1) + wni(i-1,j  ,k+1))*0.5
        wxminus = (wni(i  ,j  ,k  ) + wni(i-1,j  ,k  ))*0.5

        ! get derivatives at 1/2 locations
        dudxp = (uni(i+1,j  ,k  ) - uni(i  ,j  ,k  )) * dx(i  , CENTER    )
        dudxm = (uni(i  ,j  ,k  ) - uni(i-1,j  ,k  )) * dx(i-1, CENTER    )
        dudyp = (uni(i  ,j+1,k  ) - uni(i  ,j  ,k  )) * dy(j  , RIGHT_EDGE)
        dudym = (uni(i  ,j  ,k  ) - uni(i  ,j-1,k  )) * dy(j  , LEFT_EDGE )
        dudzp = (uni(i  ,j  ,k+1) - uni(i  ,j  ,k  )) * dz(k  , RIGHT_EDGE)
        dudzm = (uni(i  ,j  ,k  ) - uni(i  ,j  ,k-1)) * dz(k  , LEFT_EDGE )
        dvdxp = (vni(i  ,j+1,k  ) - vni(i-1,j+1,k  )) * dx(i  , LEFT_EDGE )
        dvdxm = (vni(i  ,j  ,k  ) - vni(i-1,j  ,k  )) * dx(i  , LEFT_EDGE )
        dwdxp = (wni(i  ,j  ,k+1) - wni(i-1,j  ,k+1)) * dx(i  , LEFT_EDGE )
        dwdxm = (wni(i  ,j  ,k  ) - wni(i-1,j  ,k  )) * dx(i  , LEFT_EDGE )

        ! get nu_t

        !****** requires DIAGONALS for corner ghost cells ******

        tvjp = 0.25*(tv(i-1,j  ,k  ) + tv(i  ,j  ,k  ) + &
                     tv(i  ,j+1,k  ) + tv(i-1,j+1,k  ))
        tvjm = 0.25*(tv(i-1,j-1,k  ) + tv(i  ,j-1,k  ) + &
                     tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))
        tvkp = 0.25*(tv(i-1,j  ,k  ) + tv(i  ,j  ,k  ) + &
                     tv(i  ,j  ,k+1) + tv(i-1,j  ,k+1))
        tvkm = 0.25*(tv(i-1,j  ,k-1) + tv(i  ,j,  k-1) + &
                     tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))

        ! flux of normal total stresses
        txxp = (ru1 + 2.0*tv(i,j,k))*dudxp
        txxm = (ru1 + 2.0*tv(i-1,j,k))*dudxm
        tyyp = (ru1 + tvjp)*dudyp
        tyym = (ru1 + tvjm)*dudym
        tzzp = (ru1 + tvkp)*dudzp
        tzzm = (ru1 + tvkm)*dudzm

        ! flux of cross SGS stresses
        txyp = tvjp*dvdxp
        txym = tvjm*dvdxm
        txzp = tvkp*dwdxp
        txzm = tvkm*dwdxm

        ! calculate RHS for u-momentum
        ru(i,j,k) =                                                     &                              
                   - (uxplus*uxplus - uxminus*uxminus)*dx(i,LEFT_EDGE)  &! advection term
                   - (vxplus*uyplus - vxminus*uyminus)*dy(j,CENTER)     &
                   - (wxplus*uzplus - wxminus*uzminus)*dz(k,CENTER)     &             
                   + (txxp - txxm)*dx(i,LEFT_EDGE)                      &! diffusion - normal terms
                   + (tyyp - tyym)*dy(j,CENTER)                         &
                   + (tzzp - tzzm)*dz(k,CENTER)                         &
                   + (txyp - txym)*dy(j,CENTER)                         &! diffusion - cross terms
                   + (txzp - txzm)*dz(k,CENTER)   
      enddo
    enddo
  enddo

      !++++++++++  V-COMPONENT  ++++++++++

  do k = kz1,kz2
    do j = jy1,jy2+1
      do i = ix1,ix2
 
        ! get velocities at 1/2 locations
        vxplus  = (vni(i+1,j  ,k  ) + vni(i  ,j  ,k  ))*0.5
        vxminus = (vni(i  ,j  ,k  ) + vni(i-1,j  ,k  ))*0.5

        vyplus  = (vni(i  ,j+1,k  ) + vni(i  ,j  ,k  ))*0.5
        vyminus = (vni(i  ,j  ,k  ) + vni(i  ,j-1,k  ))*0.5

        vzplus  = (vni(i  ,j  ,k+1) + vni(i  ,j  ,k  ))*0.5
        vzminus = (vni(i  ,j  ,k  ) + vni(i  ,j  ,k-1))*0.5

        uyplus  = (uni(i+1,j  ,k  ) + uni(i+1,j-1,k  ))*0.5
        uyminus = (uni(i  ,j  ,k  ) + uni(i  ,j-1,k  ))*0.5

        wyplus  = (wni(i  ,j  ,k+1) + wni(i  ,j-1,k+1))*0.5
        wyminus = (wni(i  ,j  ,k  ) + wni(i  ,j-1,k  ))*0.5

        ! get derivatives at 1/2 locations
        dvdxp = (vni(i+1,j  ,k  ) - vni(i  ,j  ,k  )) * dx(i  , RIGHT_EDGE)
        dvdxm = (vni(i  ,j  ,k  ) - vni(i-1,j  ,k  )) * dx(i  , LEFT_EDGE )
        dvdyp = (vni(i  ,j+1,k  ) - vni(i  ,j  ,k  )) * dy(j  , CENTER    )
        dvdym = (vni(i  ,j  ,k  ) - vni(i  ,j-1,k  )) * dy(j-1, CENTER    )
        dvdzp = (vni(i  ,j  ,k+1) - vni(i  ,j  ,k  )) * dz(k  , RIGHT_EDGE)
        dvdzm = (vni(i  ,j  ,k  ) - vni(i  ,j  ,k-1)) * dz(k  , LEFT_EDGE )
        dudyp = (uni(i+1,j  ,k  ) - uni(i+1,j-1,k  )) * dy(j  , LEFT_EDGE )
        dudym = (uni(i  ,j  ,k  ) - uni(i  ,j-1,k  )) * dy(j  , LEFT_EDGE )
        dwdyp = (wni(i  ,j  ,k+1) - wni(i  ,j-1,k+1)) * dy(j  , LEFT_EDGE )
        dwdym = (wni(i  ,j  ,k  ) - wni(i  ,j-1,k  )) * dy(j  , LEFT_EDGE )

        ! get nu_t

        !****** requires DIAGONALS for corner ghost cells ******

        tvip = 0.25*(tv(i  ,j-1,k  ) + tv(i+1,j-1,k  ) + &
                     tv(i+1,j  ,k  ) + tv(i  ,j  ,k  ))
        tvim = 0.25*(tv(i-1,j-1,k  ) + tv(i  ,j-1,k  ) + &
                     tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))
        tvkp = 0.25*(tv(i  ,j-1,k  ) + tv(i  ,j-1,k+1) + &
                     tv(i  ,j  ,k+1) + tv(i  ,j  ,k  ))
        tvkm = 0.25*(tv(i  ,j-1,k-1) + tv(i  ,j-1,k  ) + &
                     tv(i  ,j  ,k  ) + tv(i  ,j  ,k-1))

        ! flux of normal total stresses
        txxp = (ru1 + tvip)*dvdxp
        txxm = (ru1 + tvim)*dvdxm
        tyyp = (ru1 + 2.0*tv(i,j,k))*dvdyp
        tyym = (ru1 + 2.0*tv(i,j-1,k))*dvdym
        tzzp = (ru1 + tvkp)*dvdzp
        tzzm = (ru1 + tvkm)*dvdzm

        ! flux of cross SGS stresses
        txyp = tvip*dudyp
        txym = tvim*dudym
        tyzp = tvkp*dwdyp
        tyzm = tvkm*dwdym

        ! calculate RHS for v-momentum
        rv(i,j,k) =                                                     &
                   - (uyplus*vxplus - uyminus*vxminus)*dx(i,CENTER)     &! advection term
                   - (vyplus*vyplus - vyminus*vyminus)*dy(j,LEFT_EDGE)  &
                   - (wyplus*vzplus - wyminus*vzminus)*dz(k,CENTER)     &
                   + (txxp - txxm)*dx(i,CENTER)                         &! diffusion - normal terms
                   + (tyyp - tyym)*dy(j,LEFT_EDGE)                      &
                   + (tzzp - tzzm)*dz(k,CENTER)                         &
                   + (txyp - txym)*dx(i,CENTER)                         &! diffusion - cross terms
                   + (tyzp - tyzm)*dz(k,CENTER)  
      enddo
    enddo
  enddo

  !++++++++++  W-COMPONENT  ++++++++++
      
  do k = kz1,kz2+1
    do j = jy1,jy2
      do i = ix1,ix2

        ! get velocities at 1/2 locations
        wxplus  = (wni(i+1,j  ,k  ) + wni(i  ,j  ,k  ))*0.5
        wxminus = (wni(i  ,j  ,k  ) + wni(i-1,j  ,k  ))*0.5
               
        wyplus  = (wni(i  ,j+1,k  ) + wni(i  ,j  ,k  ))*0.5
        wyminus = (wni(i  ,j  ,k  ) + wni(i  ,j-1,k  ))*0.5
               
        wzplus  = (wni(i  ,j  ,k+1) + wni(i  ,j  ,k  ))*0.5
        wzminus = (wni(i  ,j  ,k  ) + wni(i  ,j  ,k-1))*0.5

        uzplus  = (uni(i+1,j  ,k  ) + uni(i+1,j  ,k-1))*0.5
        uzminus = (uni(i  ,j  ,k  ) + uni(i  ,j  ,k-1))*0.5

        vzplus  = (vni(i  ,j+1,k  ) + vni(i  ,j+1,k-1))*0.5
        vzminus = (vni(i  ,j  ,k  ) + vni(i  ,j  ,k-1))*0.5

        ! get temperature at z locations
        tz = (tni(i,j,k) + tni(i,j,k-1))*0.5
               
        ! get derivatives at 1/2 locations
        dwdxp = (wni(i+1,j  ,k  ) - wni(i  ,j  ,k  )) * dx(i  , RIGHT_EDGE)
        dwdxm = (wni(i  ,j  ,k  ) - wni(i-1,j  ,k  )) * dx(i  , LEFT_EDGE )
        dwdyp = (wni(i  ,j+1,k  ) - wni(i  ,j  ,k  )) * dy(j  , RIGHT_EDGE)
        dwdym = (wni(i  ,j  ,k  ) - wni(i  ,j-1,k  )) * dy(j  , LEFT_EDGE )
        dwdzp = (wni(i  ,j  ,k+1) - wni(i  ,j  ,k  )) * dz(i  , CENTER    )
        dwdzm = (wni(i  ,j  ,k  ) - wni(i  ,j  ,k-1)) * dz(i-1, CENTER    )
        dudzp = (uni(i+1,j  ,k  ) - uni(i+1,j  ,k-1)) * dz(k  , LEFT_EDGE )
        dudzm = (uni(i  ,j  ,k  ) - uni(i  ,j  ,k-1)) * dz(k  , LEFT_EDGE )
        dvdzp = (vni(i  ,j+1,k  ) - vni(i  ,j+1,k-1)) * dz(k  , LEFT_EDGE )
        dvdzm = (vni(i  ,j  ,k  ) - vni(i  ,j  ,k-1)) * dz(k  , LEFT_EDGE )

        ! get nu_t

        !****** requires DIAGONALS for corner ghost cells ******

        tvip = 0.25*(tv(i  ,j  ,k-1) + tv(i+1,j  ,k-1) + &
                     tv(i+1,j  ,k  ) + tv(i  ,j  ,k  ))
        tvim = 0.25*(tv(i-1,j  ,k-1) + tv(i  ,j  ,k-1) + &
                     tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))
        tvjp = 0.25*(tv(i  ,j  ,k-1) + tv(i  ,j  ,k  ) + &
                     tv(i  ,j+1,k  ) + tv(i  ,j+1,k-1))
        tvjm = 0.25*(tv(i  ,j-1,k-1) + tv(i  ,j-1,k  ) + & 
                     tv(i  ,j  ,k  ) + tv(i  ,j  ,k-1))

        ! flux of normal total stresses
        txxp = (ru1 + tvip)*dwdxp
        txxm = (ru1 + tvim)*dwdxm
        tyyp = (ru1 + tvjp)*dwdyp
        tyym = (ru1 + tvjm)*dwdym
        tzzp = (ru1 + 2.0*tv(i,j,k))*dwdzp
        tzzm = (ru1 + 2.0*tv(i,j,k-1))*dwdzm

        ! flux of cross SGS stresses
        txzp = tvip*dudzp
        txzm = tvim*dudzm
        tyzp = tvjp*dvdzp
        tyzm = tvjm*dvdzm

        ! calculate RHS for w-momentum
        rw(i,j,k) =                                                     &
                   - (uzplus*wxplus - uzminus*wxminus)*dx(i,CENTER)     &! advection term
                   - (vzplus*wyplus - vzminus*wyminus)*dy(j,CENTER)     &
                   - (wzplus*wzplus - wzminus*wzminus)*dz(k,LEFT_EDGE)  &
                   + (txxp - txxm)*dx(i,CENTER)                         &! diffusion - normal terms
                   + (tyyp - tyym)*dy(j,CENTER)                         &
                   + (tzzp - tzzm)*dz(k,LEFT_EDGE)                      &
                   + (txzp - txzm)*dx(i,CENTER)                         &! diffusion - cross terms
                   + (tyzp - tyzm)*dy(j,CENTER)                         & 
                   + tz * tz_m                                           ! energy momentum coupling
      enddo
    enddo
  enddo

  END SUBROUTINE ins_rhs3d

