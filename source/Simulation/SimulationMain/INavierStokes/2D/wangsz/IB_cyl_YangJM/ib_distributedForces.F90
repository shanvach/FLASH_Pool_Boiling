

#include "constants.h"
#include "Flash.h"
#include "ImBound.h"

!#define TWO_POINTSP 1
!#define TWO_NORMAL_POINTS 1
#define TWO_NORMAL_POINTS_EUL 1
!#define PRES_GRAD_CORR 1
!#define PRES_BDL_EQ 1
!#define PRES_NS_EQ  12

subroutine ib_distributedForces(blockID, particleData, vortx, vorty, vortz)

  use Grid_Data, only : gr_meshMe

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,      &
                             Grid_getDeltas, Grid_getBlkCenterCoords, &
                             Grid_getBlkPhysicalSize

  use ImBound_Data, only : ib_nu,ib_stencil,ib_alphax,ib_alphay,ib_alphaz,ib_dt  !, &
!                           ib_alpha

  use ib_interface, only : ib_stencils,ib_getInterpFunc

  use Driver_interface, only : Driver_abortFlash

  use Grid_data,ONLY : gr_imin,gr_jmin,gr_kmin

  use IncompNS_data, ONLY : ins_gravX,ins_gravY,ins_gravZ

  implicit none
  integer, intent(IN) :: blockID
  real, intent(IN), dimension(GRID_IHI_GC*K1D+1,GRID_JHI_GC*K2D+1,GRID_KHI_GC*K3D+1) :: vortx 
  real, intent(IN), dimension(GRID_IHI_GC*K1D+1,GRID_JHI_GC*K2D+1,GRID_KHI_GC*K3D+1) :: vorty
  real, intent(IN), dimension(GRID_IHI_GC*K1D+1,GRID_JHI_GC*K2D+1,GRID_KHI_GC*K3D+1) :: vortz
  real, intent(INOUT) :: particleData(NPART_PROPS)

  ! Local Variables....
  real :: xp,yp,zp,zL,h,hl,dx,dy,dz,dsx,dsy,dsz,ubd,vbd,wbd,ubdd,vbdd,wbdd,nxp,nyp,nzp
  real, dimension(MDIM) :: xbe,del,coord,bsize,np

  integer, parameter, dimension(MDIM)      :: grdip = (/ CENTER, CENTER, CENTER /)
  real, parameter, dimension(MDIM)      :: dlip = (/ 0.5, 0.5, 0.5 /)


#ifdef TANGENT_WITH_VORTICITY
  integer, parameter, dimension(MDIM,MDIM) :: grdnu = &
           RESHAPE( (/  FACES, FACES,CENTER, CENTER, FACES, FACES, FACES,CENTER, FACES /), (/MDIM,MDIM /))
                    !    k (wz)               i (wx)                j (wy)
  real, parameter, dimension(MDIM,MDIM) :: dlnu = &
           RESHAPE( (/ 0., 0., 0.5, 0.5, 0., 0., 0., 0.5, 0. /), (/MDIM,MDIM /))  
                    !   k             i           j
#else
  integer, parameter, dimension(MDIM,MDIM) :: grdu = &
           RESHAPE( (/  FACES,CENTER,CENTER, CENTER, FACES,CENTER, CENTER,CENTER, FACES /), (/MDIM,MDIM /))
                    !    u                   v                     w
  real, parameter, dimension(MDIM,MDIM) :: dlu = &
           RESHAPE( (/ 0., 0.5, 0.5, 0.5, 0., 0.5, 0.5, 0.5, 0. /), (/MDIM,MDIM /))
                    !  u             v             w
  real :: ui,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,exx,eyy,ezz,exy,exz,eyz
  real :: ue,ve,we,ues,ves,wes,ven,ups,vps,wps,vpn,dun,dvn,dwn,dune,dvne,dwne

  real :: normt,tx,ty,tz,dpdxt,dpdx,dpdy,dpdz

  real :: zpres2,xbe2(MDIM),zpres3,xbe3(MDIM)

  ! Shizhao Wang for testing the pressure correction
  real :: dudx2,dudy2,dudz2,dvdx2,dvdy2,dvdz2,dwdx2,dwdy2,dwdz2,exx2,eyy2,ezz2,exy2,exz2,eyz2
  real :: ue2,ve2,we2,ues2,ves2,wes2,ven2
  real :: dutdxn_e1, dutdxn_e2, dun_e, dvn_e, dwn_e
  real :: ddun_e, ddvn_e, ddwn_e
  real :: dun_2n_l, dvn_2n_l, dun_2n_vg, dvn_2n_vg
  real :: zL2
  real :: wz_2n_l, wz_2n_vg, wz_2n_e, wz_2n_e1, wz_2n_e2
  real :: dpdxt2,dpdx2,dpdy2,dpdz2, dpdn2
  real :: wz_2n_lw, wz_2n_le

  real :: pppx, pppy, pppxpx, pppxpy, pppypx, pppypy
  real :: pupx, pupy, pvpx, pvpy
  real :: pupxpx, pupxpy, pupypx, pupypy, pvpxpx, pvpxpy, pvpypx, pvpypy
  real :: dpdxtdxn
  real :: dun_l_vg, dvn_l_vg, dun_l_2p, dvn_l_2p
  real :: dun_c, dvn_c, dun_c2, dvn_c2
  real :: dun_vgc, dvn_vgc, dun_vgc2, dvn_vgc2
  real :: dun_2pc, dvn_2pc, dun_2pc2, dvn_2pc2
  real :: wz_l_vg, wz_l_2p, wz_c, wz_c2
  real :: wz_vgc, wz_vgc2, wz_2pc, wz_2pc2
 
#endif

  integer :: presflag,gridfl(MDIM),i,idim,gridind,nkij
  integer :: ielem(ib_stencil,MDIM,CONSTANT_TWO)
  real :: dpdn,eps
  real :: delaux(MDIM),xyz_stencil(ib_stencil,MDIM),phile(ib_stencil,NDIM+1)

  real :: zpres,p_i,zv(MDIM),nuwx,nuwy,nuwz

  real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData,facezData

#ifdef TWO_NORMAL_POINTS_EUL
  real :: enh1, enh2, enh3
  real :: fvxen1, fvxen2, fvxw
  real :: fvyen1, fvyen2, fvyw
  real :: sxxen1, sxyen1, syxen1, syyen1
  real :: sxxen2, sxyen2, syxen2, syyen2
  real :: uen1, uen1xp, uen1xm, uen1yp, uen1ym
  real :: ven1, ven1xp, ven1xm, ven1yp, ven1ym
  real :: uen2, uen2xp, uen2xm, uen2yp, uen2ym
  real ::       ven2xp, ven2xm, ven2yp, ven2ym
  real :: uixp, uixm,uiyp, uiym
  real :: uixpyp, uixpym, uixmyp, uixmym
  real :: ddudxdx, ddudxdy, ddudydx, ddudydy, ddvdxdx, ddvdxdy, ddvdydx, ddvdydy
  real :: ddudxdx2, ddudxdy2, ddudydx2, ddudydy2, ddvdxdx2, ddvdxdy2, ddvdydx2, ddvdydy2
  real :: fvxen1_c, fvyen1_c, fvxw_c, fvyw_c, wz_2n_e1_c, wz_2n_w_c
  real :: fvxen2_c, fvyen2_c, fvxw_c2, fvyw_c2, wz_2n_e2_c, wz_2n_w_c2
  real :: pres_c2, pppxn2
  real :: pixp, pixm, piyp, piym
  real :: pppxt2, ppppxnpxt2,  ppppxpx, ppppxpy, ppppypy
  real :: pixpyp, pixpym, pixmyp, pixmym
  real :: wz_2n_e2_pc, wz_2n_e2_ppc
  real :: fvxen2_pc, fvyen2_pc, fvxen2_ppc, fvyen2_ppc
  real :: fvxen2_c2, fvyen2_c2, wz_2n_e2_c2
  real :: pppxn1, pppxt1
  real :: ppppxnpxt3,  ppppxpx3, ppppxpy3, ppppypy3
  real :: wz_2n_e2_ppc3
  real :: fvxen2_ppc3, fvyen2_ppc3
 
#endif

#ifdef TANGENT_WITH_VORTICITY
  integer, parameter :: derivflag = 0 ! Give interpolation functions.
#else
  integer, parameter :: derivflag = 1 ! Give interpolation functions and their derivatives
#endif

  real ::nu, dt
  real :: pi, x, y

!  write(*,*) 'ib_alpha', ib_alpha

  pi = acos(-1.0)
!  ib_alphax = ib_alpha
!  ib_alphay = ib_alpha
!  ib_alphax = ib_alpha
  
  nu = ib_nu
  dt = ib_dt

#ifndef TANGENT_WITH_VORTICITY  ! Shizhao Wang
  ue =0.; ve=0.; we=0.; 
  ue2 =0.; ve2=0.; we2=0.; 

  dudx=0.; dudy=0.; dudz=0.;
  dvdx=0.; dvdy=0.; dvdz=0.; 
  dwdx=0.; dwdy=0.; dwdz=0.;

  dudx2=0.; dudy2=0.; dudz2=0.;
  dvdx2=0.; dvdy2=0.; dvdz2=0.; 
  dwdx2=0.; dwdy2=0.; dwdz2=0.;

  exx=0;  exy=0.; exz=0.;
  eyy=0.; eyz=0.;
  ezz=0.;
#endif

  ! Get dx,dy
  call Grid_getDeltas(blockID,del)
  call Grid_getBlkCenterCoords(blockID,coord)
  call Grid_getBlkPhysicalSize(blockID,bsize)

  dx=Del(IAXIS)
  dy=Del(JAXIS)
#if NDIM == 3
  dz=Del(KAXIS)
#else
  dz=1.
#endif

  ! Particle data:
  xp  = particleData(POSX_PART_PROP)
  yp  = particleData(POSY_PART_PROP)
  ubd = particleData(VELX_PART_PROP)
  vbd = particleData(VELY_PART_PROP)
  ubdd= particleData(ACCX_PART_PROP)
  vbdd= particleData(ACCY_PART_PROP)
  nxp = particleData(NMLX_PART_PROP)
  nyp = particleData(NMLY_PART_PROP)

  np(IAXIS) = nxp
  np(JAXIS) = nyp

#if NDIM == 3
  zp  = particleData(POSZ_PART_PROP)
  wbd = particleData(VELZ_PART_PROP)
  wbdd= particleData(ACCZ_PART_PROP)
  nzp = particleData(NMLZ_PART_PROP)
#else
  zp  = 0.
  wbd = 0.
  wbdd= 0.
  nzp = 0.
#endif
  np(KAXIS) = nzp


  ! Function Gimme h:
  ! call ib_normaldistance(dx,dy,dz,h)
  dsx = ib_alphax*dx
  dsy = ib_alphay*dy
#if NDIM == 2
  !h   = 0.5*(dsx+dsy) 
  h   = 1.*sqrt( (dsx*nxp)**2. + (dsy*nyp)**2. )
  eps = 1.E-2*MIN(dx,dy)  
#elif NDIM == 3
  dsz = ib_alphaz*dz
  !h   = 1./3.*(dsx+dsy+dsz)
  h   = 1.*sqrt( (dsx*nxp)**2. + (dsy*nyp)**2. + (dsz*nzp)**2. )
  eps = 1.E-2*MIN(dx,dy,dz)
#endif

             
  ! External Point Position:
  xbe(IAXIS) = xp + nxp*h
  xbe(JAXIS) = yp + nyp*h
#if NDIM == 3
  xbe(KAXIS) = zp + nzp*h
#else
  xbe(KAXIS) = 0.
#endif

  zpres = 0.
  zv(1:MDIM) = 0.
  zL= 0.
  nuwx = 0.
  nuwy = 0.
  nuwz = 0.
  do presflag = CONSTANT_ZERO,CONSTANT_ONE

     ! N kij
#ifdef TANGENT_WITH_VORTICITY
     nkij = 1 + (1-presflag)*(2*NDIM-MDIM-1)
#else
     nkij = 1 + (1-presflag)*(NDIM-1) ! in 2D d/dx,d/dy; in 3D d/dx,d/dy,d/dz
#endif

     do gridind = 1,nkij

#ifdef TANGENT_WITH_VORTICITY
        ! Define Grids in case of vorticity and pressure:
        gridfl(:) = presflag*grdip(:) + (1-presflag)*grdnu(:,gridind)
        ! Auxiliary deltas
        do idim = 1,NDIM
           delaux(idim) = (real(presflag)*dlip(idim)+real(1-presflag)*dlnu(idim,gridind))*del(idim)
        enddo
#else
        ! Define Grids in case of vorticity and pressure:
        gridfl(:) = presflag*grdip(:) + (1-presflag)*grdu(:,gridind)
        ! Auxiliary deltas
        do idim = 1,NDIM
           delaux(idim) = (real(presflag)*dlip(idim)+real(1-presflag)*dlu(idim,gridind))*del(idim)
        enddo
#endif

        ! Obtain Stencil for External Point:
        call ib_stencils(xbe,np,gridfl,del,coord,bsize,   & 
                         ielem(:,:,presflag+1),hl,COMPUTE_FORCES)


        ! Compute shape functions
        ! Positions of points on the stencil:
        xyz_stencil(1:ib_stencil,1:MDIM) = 0. 
        do idim = 1,NDIM
           xyz_stencil(1:ib_stencil,idim) = coord(idim) - 0.5*bsize(idim) + &
                real(ielem(1:ib_stencil,idim,presflag+1) - NGUARD - 1)*del(idim) + delaux(idim) 
        enddo

        ! Get interpolation functions:
        call ib_getInterpFunc(xbe,xyz_stencil,del,derivflag,phile)

!  print*, 'ib_stencil, gridind, presflag', ib_stencil, gridind, presflag
!  print*, 'stencil', phile

        if (presflag .eq. CONSTANT_ONE) then         ! Pressure

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,solnData,CENTER)

#ifndef TANGENT_WITH_VORTICITY
           dpdx = 0.; dpdy =0.;
           dpdz = 0.
#endif
           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              p_i = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                      ielem(i,JAXIS,presflag+1), &
                                      ielem(i,KAXIS,presflag+1));  

              zpres = zpres + phile(i,1)*p_i 

#ifndef TANGENT_WITH_VORTICITY
              dpdx = dpdx + phile(i,2)*p_i;
              dpdy = dpdy + phile(i,3)*p_i;
#if NDIM == MDIM
              dpdz = dpdz + phile(i,4)*p_i;
#endif             
#endif
           enddo

           ! Release Pointer
           call Grid_releaseBlkPtr(blockID,solnData,CENTER)

           ! Get Pressure approximation at surface marker: Acceleration +
           ! gravity effects.
           dpdn = -(     ubdd*nxp +      vbdd*nyp +      wbdd*nzp) + & ! -rho*Du/Dt * n 
                   (ins_gravX*nxp + ins_gravY*nyp + ins_gravZ*nzp);    ! +rho*    g * n
           zL = zpres - dpdn*h;

        else                                         ! Tangent stress

#ifdef TANGENT_WITH_VORTICITY
 
           select case (gridind)
           case(1) 
           ! Only wz:
           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              zv(gridind) = zv(gridind) + phile(i,1)*vortz(ielem(i,IAXIS,presflag+1), &
                                                           ielem(i,JAXIS,presflag+1), &
                                                           ielem(i,KAXIS,presflag+1)); 
           enddo
           nuwz = nu*zv(gridind)

           case(2)
           ! wx:
           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              zv(gridind) = zv(gridind) + phile(i,1)*vortx(ielem(i,IAXIS,presflag+1), &
                                                           ielem(i,JAXIS,presflag+1), &
                                                           ielem(i,KAXIS,presflag+1)); 
           enddo
           nuwx = nu*zv(gridind)          

           case(3)
           ! wy:
           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              zv(gridind) = zv(gridind) + phile(i,1)*vorty(ielem(i,IAXIS,presflag+1), &
                                                           ielem(i,JAXIS,presflag+1), &
                                                           ielem(i,KAXIS,presflag+1)); 
           enddo
           nuwy = nu*zv(gridind)

           end select

#else

           select case(gridind)
           case(1)

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,facexData,FACEX)

           ! U velocity derivatives on external point:
           ue = 0.; dudx = 0.; dudy =0.;
           dudz = 0.
           do i = 1 , ib_stencil

              ui = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              ue   = ue   + phile(i,1)*ui;
              dudx = dudx + phile(i,2)*ui;
              dudy = dudy + phile(i,3)*ui;
#if NDIM == MDIM
              dudz = dudz + phile(i,4)*ui;
#endif             
           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,facexData,FACEX)

           case(2)

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,faceyData,FACEY)

           ! V velocity derivatives on external point:
           ve = 0.; dvdx = 0.; dvdy =0.;
           dvdz = 0.
           do i = 1 , ib_stencil

              ui = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              ve   = ve   + phile(i,1)*ui; 
              dvdx = dvdx + phile(i,2)*ui;
              dvdy = dvdy + phile(i,3)*ui;
#if NDIM == MDIM
              dvdz = dvdz + phile(i,4)*ui;
#endif             
           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,faceyData,FACEY)


           case(3) 

#if NDIM == MDIM
           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,facezData,FACEZ)

           ! W velocity derivatives on external point:
           we   = 0.; dwdx = 0.; dwdy =0.;
           dwdz = 0.
           do i = 1 , ib_stencil

              ui = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              we   = we   + phile(i,1)*ui;
              dwdx = dwdx + phile(i,2)*ui;
              dwdy = dwdy + phile(i,3)*ui;
              dWdz = dwdz + phile(i,4)*ui;

           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

           end select

#endif
           
        end if

     enddo

  enddo

!######################################################################################################################

#ifdef TWO_NORMAL_POINTS
  ! External Point Position:
  xbe2(IAXIS) = xp + nxp*h*2
  xbe2(JAXIS) = yp + nyp*h*2
#if NDIM == 3
  xbe2(KAXIS) = zp + nzp*h*2
#else
  xbe2(KAXIS) = 0.
#endif

  zpres2 = 0.
!  zv2(1:MDIM) = 0.
  zL2= 0.
  do presflag = CONSTANT_ZERO,CONSTANT_ONE

     ! N kij
#ifdef TANGENT_WITH_VORTICITY
     nkij = 1 + (1-presflag)*(2*NDIM-MDIM-1)
#else
     nkij = 1 + (1-presflag)*(NDIM-1) ! in 2D d/dx,d/dy; in 3D d/dx,d/dy,d/dz
#endif

     do gridind = 1,nkij

#ifdef TANGENT_WITH_VORTICITY
        print*, 'TWO_NORMAL_POINTS does not work for TANGENT_WITH_VORTICITY CASE'
        pause
#else
        ! Define Grids in case of vorticity and pressure:
        gridfl(:) = presflag*grdip(:) + (1-presflag)*grdu(:,gridind)
        ! Auxiliary deltas
        do idim = 1,NDIM
           delaux(idim) = (real(presflag)*dlip(idim)+real(1-presflag)*dlu(idim,gridind))*del(idim)
        enddo
#endif

        ! Obtain Stencil for External Point:
        call ib_stencils(xbe2,np,gridfl,del,coord,bsize,   & 
                         ielem(:,:,presflag+1),hl,COMPUTE_FORCES)


        ! Compute shape functions
        ! Positions of points on the stencil:
        xyz_stencil(1:ib_stencil,1:MDIM) = 0. 
        do idim = 1,NDIM
           xyz_stencil(1:ib_stencil,idim) = coord(idim) - 0.5*bsize(idim) + &
                real(ielem(1:ib_stencil,idim,presflag+1) - NGUARD - 1)*del(idim) + delaux(idim) 
        enddo

        ! Get interpolation functions:
        call ib_getInterpFunc(xbe2,xyz_stencil,del,derivflag,phile)

!  print*, 'ib_stencil, gridind, presflag', ib_stencil, gridind, presflag
!  print*, 'stencil', phile

        if (presflag .eq. CONSTANT_ONE) then         ! Pressure

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,solnData,CENTER)

#ifndef TANGENT_WITH_VORTICITY
           dpdx2 = 0.; dpdy2 =0.;
           dpdz2 = 0.
#endif
           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              p_i = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                      ielem(i,JAXIS,presflag+1), &
                                      ielem(i,KAXIS,presflag+1));  

              zpres2 = zpres2 + phile(i,1)*p_i 

#ifndef TANGENT_WITH_VORTICITY
              dpdx2 = dpdx2 + phile(i,2)*p_i;
              dpdy2 = dpdy2 + phile(i,3)*p_i;
#if NDIM == MDIM
              dpdz2 = dpdz2 + phile(i,4)*p_i;
#endif             
#endif
           enddo

           ! Release Pointer
           call Grid_releaseBlkPtr(blockID,solnData,CENTER)

           ! Get Pressure approximation at surface marker: Acceleration +
           ! gravity effects.
           dpdn2 = -(     ubdd*nxp +      vbdd*nyp +      wbdd*nzp) + & ! -rho*Du/Dt * n 
                   (ins_gravX*nxp + ins_gravY*nyp + ins_gravZ*nzp);    ! +rho*    g * n
           zL2 = zpres2 - dpdn*h*2;

        else                                         ! Tangent stress

#ifdef TANGENT_WITH_VORTICITY
 
        print*, 'TWO_NORMAL_POINTS does not work for TANGENT_WITH_VORTICITY CASE'
#else

           select case(gridind)
           case(1)

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,facexData,FACEX)

           ! U velocity derivatives on external point:
           ue2 = 0.; dudx2 = 0.; dudy2 =0.;
           dudz2 = 0.
           do i = 1 , ib_stencil

              ui = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              ue2   = ue2   + phile(i,1)*ui;
              dudx2 = dudx2 + phile(i,2)*ui;
              dudy2 = dudy2 + phile(i,3)*ui;
#if NDIM == MDIM
              dudz2 = dudz2 + phile(i,4)*ui;
#endif             
           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,facexData,FACEX)

           case(2)

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,faceyData,FACEY)

           ! V velocity derivatives on external point:
           ve2 = 0.; dvdx2 = 0.; dvdy2 =0.;
           dvdz2 = 0.
           do i = 1 , ib_stencil

              ui = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              ve2   = ve2   + phile(i,1)*ui; 
              dvdx2 = dvdx2 + phile(i,2)*ui;
              dvdy2 = dvdy2 + phile(i,3)*ui;
#if NDIM == MDIM
              dvdz2 = dvdz2 + phile(i,4)*ui;
#endif             
           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,faceyData,FACEY)


           case(3) 

#if NDIM == MDIM
           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,facezData,FACEZ)

           ! W velocity derivatives on external point:
           we2   = 0.; dwdx2 = 0.; dwdy2 =0.;
           dwdz2 = 0.
           do i = 1 , ib_stencil

              ui = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              we2   = we2   + phile(i,1)*ui;
              dwdx2 = dwdx2 + phile(i,2)*ui;
              dwdy2 = dwdy2 + phile(i,3)*ui;
              dWdz2 = dwdz2 + phile(i,4)*ui;

           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

           end select

#endif
           
        end if

     enddo

  enddo

  ! Project external and marker Velocity on the plane:
  ven = nxp*ue + nyp*ve + nzp*we

  ! ves = ve - (ve*n) n
  ues  = ue - ven*nxp
  ves  = ve - ven*nyp
  wes  = we - ven*nzp

  ! Project external and marker Velocity on the plane:
  ven2 = nxp*ue2 + nyp*ve2 + nzp*we2

  ! ves = ve - (ve*n) n
  ues2  = ue2 - ven2*nxp
  ves2  = ve2 - ven2*nyp
  wes2  = we2 - ven2*nzp

  ! Linear Part:
  normt = 1./h
  dun_e = (ues2-ues)*normt
  dvn_e = (ves2-ves)*normt
  dwn_e = (wes2-wes)*normt

  normt = 1./sqrt(dun_e**2. + dvn_e**2. + dwn_e**2.)
  tx = dun_e*normt
  ty = dvn_e*normt
  tz = dwn_e*normt

  ! slope
  if (NDIM==MDIM) pause 'slope not defined for 3D flows'
  dutdxn_e1 = nxp*dudx*tx + nxp*dvdx*ty + nyp*dudy*tx + nyp*dvdy*ty 
  dutdxn_e2 = nxp*dudx2*tx + nxp*dvdx2*ty + nyp*dudy2*tx + nyp*dvdy2*ty 

  ddun_e = (dutdxn_e2-dutdxn_e1)/h*tx
  ddvn_e = (dutdxn_e2-dutdxn_e1)/h*ty

  dun_2n_l = dun_e - ddun_e*h
  dvn_2n_l = dvn_e - ddvn_e*h

  dun_2n_vg = (2.*dutdxn_e1-dutdxn_e2)*tx
  dvn_2n_vg = (2.*dutdxn_e1-dutdxn_e2)*ty

  wz_2n_l  = dvn_2n_l*nxp - dun_2n_l*nyp
  wz_2n_vg = dvn_2n_vg*nxp - dun_2n_vg*nyp
  
  wz_2n_e = dvn_e*nxp - dun_e*nyp
  wz_2n_e1 = dutdxn_e1*ty*nxp - dutdxn_e1*tx*nyp
  wz_2n_e2 = dutdxn_e2*ty*nxp - dutdxn_e2*tx*nyp
   
  write(9500+gr_meshMe, '(30f20.8)') particleData(GLOB_PART_PROP), zL, zL2,  &
                       nu*dun_2n_l, nu*dvn_2n_l, wz_2n_l,                    &
                       nu*dun_2n_vg, nu*dvn_2n_vg, wz_2n_vg, dun_e, dvn_e,   &
                       dutdxn_e1, dutdxn_e2, wz_2n_e, wz_2n_e1, wz_2n_e2


#endif /* TWO_NORMAL_POINTS */

!######################################################################################################################

#ifdef TWO_NORMAL_POINTS_EUL

!#ifdef TWO_NORMAL_POINTS_EUL
!  real :: enh1, enh2
!  real :: fven1, fven2, fvw
!  real :: sxxen1, sxyen1, syxen1, syyen1
!  real :: sxxen2, sxyen2, syxen2, syyen2
!  real :: uen1, uen1xp, uen1xm, uen1yp, uen1ym
!  real :: ven1, ven1xp, ven1xm, ven1yp, ven1ym
!  real :: uen2, uen2xp, uen2xm, uen2yp, uen2ym
!  real :: ven2, ven2xp, ven2xm, ven2yp, ven2ym
!#endif

  enh1 = h
  !enh2 = 1.67*h
  enh2 = 2.0*dx
  !enh2 = 2.0*h
  enh3 = h/1.2*4.0
 
  fvxen1 = 0.
  fvxen2 = 0. 
  fvxw   = 0.
  fvyen1 = 0.
  fvyen2 = 0. 
  fvyw   = 0.
  sxxen1= 0. 
  sxyen1= 0. 
  syxen1= 0. 
  syyen1= 0.

  sxxen2=0.
  sxyen2=0.
  syxen2=0.
  syyen2=0.

  uen1 = 0.
  uen1xp = 0.
  uen1xm = 0. 
  uen1yp = 0. 
  uen1ym = 0.

  ven1 = 0.
  ven1xp = 0.
  ven1xm = 0. 
  ven1yp = 0. 
  ven1ym = 0.

  uen2 = 0.
  uen2xp = 0.
  uen2xm = 0.
  uen2yp = 0.
  uen2ym = 0.

  ven2 = 0.
  ven2xp = 0.
  ven2xm = 0.
  ven2yp = 0.
  ven2ym = 0.

  ! first point
  ! External Point Position:
  xbe(IAXIS) = xp + nxp*enh1
  xbe(JAXIS) = yp + nyp*enh1
#if NDIM == 3
  xbe(KAXIS) = zp + nzp*enh1
#else
  xbe(KAXIS) = 0.
#endif

  zpres = 0.
  zL= 0.
  do presflag = CONSTANT_ZERO,CONSTANT_ONE

     ! N kij
#ifdef TANGENT_WITH_VORTICITY
        print*, 'TWO_NORMAL_POINTS does not work for TANGENT_WITH_VORTICITY CASE'
#else
     nkij = 1 + (1-presflag)*(NDIM-1) ! in 2D d/dx,d/dy; in 3D d/dx,d/dy,d/dz
#endif

     do gridind = 1,nkij

#ifdef TANGENT_WITH_VORTICITY
        print*, 'TWO_NORMAL_POINTS does not work for TANGENT_WITH_VORTICITY CASE'
        pause
#else
        ! Define Grids in case of vorticity and pressure:
        gridfl(:) = presflag*grdip(:) + (1-presflag)*grdu(:,gridind)
        ! Auxiliary deltas
        do idim = 1,NDIM
           delaux(idim) = (real(presflag)*dlip(idim)+real(1-presflag)*dlu(idim,gridind))*del(idim)
        enddo
#endif

        ! Obtain Stencil for External Point:
        call ib_stencils(xbe,np,gridfl,del,coord,bsize,   & 
                         ielem(:,:,presflag+1),hl,COMPUTE_FORCES)

        ! Compute shape functions
        ! Positions of points on the stencil:
        xyz_stencil(1:ib_stencil,1:MDIM) = 0. 
        do idim = 1,NDIM
           xyz_stencil(1:ib_stencil,idim) = coord(idim) - 0.5*bsize(idim) + &
                real(ielem(1:ib_stencil,idim,presflag+1) - NGUARD - 1)*del(idim) + delaux(idim) 
        enddo

        ! Get interpolation functions:
        call ib_getInterpFunc(xbe,xyz_stencil,del,derivflag,phile)

        if (presflag .eq. CONSTANT_ONE) then         ! Pressure

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,solnData,CENTER)

#ifndef TANGENT_WITH_VORTICITY
           dpdx = 0.; dpdy =0.;
           dpdz = 0.
#endif
           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              p_i = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                      ielem(i,JAXIS,presflag+1), &
                                      ielem(i,KAXIS,presflag+1));  

              pixp = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                       ielem(i,JAXIS,presflag+1), &
                                       ielem(i,KAXIS,presflag+1));
              pixm = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                       ielem(i,JAXIS,presflag+1), &
                                       ielem(i,KAXIS,presflag+1));
              piyp = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                       ielem(i,JAXIS,presflag+1)+1, &
                                       ielem(i,KAXIS,presflag+1));
              piym = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                       ielem(i,JAXIS,presflag+1)-1, &
                                       ielem(i,KAXIS,presflag+1));

              zpres = zpres + phile(i,1)*p_i 
              dpdx  = dpdx  + phile(i,1)*(pixp-pixm)/(2.*dx)
              dpdy  = dpdy  + phile(i,1)*(piyp-piym)/(2.*dy)
#ifndef TANGENT_WITH_VORTICITY
!              dpdx = dpdx + phile(i,2)*p_i;
!              dpdy = dpdy + phile(i,3)*p_i;
#if NDIM == MDIM
              dpdz = dpdz + phile(i,4)*p_i;
#endif             
#endif
           enddo

           ! Release Pointer
           call Grid_releaseBlkPtr(blockID,solnData,CENTER)

           ! Get Pressure approximation at surface marker: Acceleration +
           ! gravity effects.
           dpdn = -(     ubdd*nxp +      vbdd*nyp +      wbdd*nzp) + & ! -rho*Du/Dt * n 
                   (ins_gravX*nxp + ins_gravY*nyp + ins_gravZ*nzp);    ! +rho*    g * n
           zL = zpres - dpdn*enh1;

        else                                         ! Tangent stress

#ifdef TANGENT_WITH_VORTICITY
 
        print*, 'TWO_NORMAL_POINTS does not work for TANGENT_WITH_VORTICITY CASE'
#else

           select case(gridind)
           case(1)

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,facexData,FACEX)

           ! U velocity derivatives on external point:
           uen1 = 0.; dudx = 0.; dudy =0.;
           dudz = 0.
           ddudxdx = 0.
           ddudxdy = 0.
           ddudydy = 0.
           do i = 1 , ib_stencil

              ui = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              uixp = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              uixm = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              uiyp = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1)+1, &
                                           ielem(i,KAXIS,presflag+1));

              uiym = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1)-1, &
                                           ielem(i,KAXIS,presflag+1));

            uixpyp = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                           ielem(i,JAXIS,presflag+1)+1, &
                                           ielem(i,KAXIS,presflag+1));

            uixmyp = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                           ielem(i,JAXIS,presflag+1)+1, &
                                           ielem(i,KAXIS,presflag+1));

            uixpym = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                           ielem(i,JAXIS,presflag+1)-1, &
                                           ielem(i,KAXIS,presflag+1));

            uixmym = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                           ielem(i,JAXIS,presflag+1)-1, &
                                           ielem(i,KAXIS,presflag+1));

              uen1   = uen1   + phile(i,1)*ui
              dudx  = dudx  + phile(i,1)*(uixp-uixm)/(2.*dx)
              dudy  = dudy  + phile(i,1)*(uiyp-uiym)/(2.*dy)
              ddudxdx = ddudxdx + phile(i,1)*(uixp+uixm-2.*ui)/(dx*dx)
              ddudydy = ddudydy + phile(i,1)*(uiyp+uiym-2.*ui)/(dy*dy)
              ddudxdy = ddudxdy + & 
                      & phile(i,1)*(uixpyp-uixpym-uixmyp+uixmym)/(4.*dx*dy)
!              ue2   = ue2   + phile(i,1)*ui;
!              dudx2 = dudx2 + phile(i,2)*ui;
!              dudy2 = dudy2 + phile(i,3)*ui;
#if NDIM == MDIM
              dudz = dudz + phile(i,4)*ui;
#endif             
           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,facexData,FACEX)

           case(2)

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,faceyData,FACEY)

           ! V velocity derivatives on external point:
           ven1 = 0.; dvdx = 0.; dvdy =0.;
           dvdz = 0.
           ddvdxdx = 0.
           ddvdxdy = 0.
           ddvdydy = 0.
           do i = 1 , ib_stencil

              ui = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              uixp = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              uixm = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              uiyp = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1)+1, &
                                           ielem(i,KAXIS,presflag+1));

              uiym = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1)-1, &
                                           ielem(i,KAXIS,presflag+1));

            uixpyp = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                           ielem(i,JAXIS,presflag+1)+1, &
                                           ielem(i,KAXIS,presflag+1));

            uixmyp = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                           ielem(i,JAXIS,presflag+1)+1, &
                                           ielem(i,KAXIS,presflag+1));

            uixpym = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                           ielem(i,JAXIS,presflag+1)-1, &
                                           ielem(i,KAXIS,presflag+1));

            uixmym = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                           ielem(i,JAXIS,presflag+1)-1, &
                                           ielem(i,KAXIS,presflag+1));

              ven1   = ven1   + phile(i,1)*ui
              dvdx  = dvdx  + phile(i,1)*(uixp-uixm)/(2.*dx)
              dvdy  = dvdy  + phile(i,1)*(uiyp-uiym)/(2.*dy)
              ddvdxdx = ddvdxdx + phile(i,1)*(uixp+uixm-2.*ui)/(dx*dx)
              ddvdydy = ddvdydy + phile(i,1)*(uiyp+uiym-2.*ui)/(dy*dy)
              ddvdxdy = ddvdxdy + & 
                      & phile(i,1)*(uixpyp-uixpym-uixmyp+uixmym)/(4.*dx*dy)
!              ve2   = ve2   + phile(i,1)*ui; 
!              dvdx2 = dvdx2 + phile(i,2)*ui;
!              dvdy2 = dvdy2 + phile(i,3)*ui;
#if NDIM == MDIM
              dvdz = dvdz + phile(i,4)*ui;
#endif             
           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,faceyData,FACEY)


           case(3) 

#if NDIM == MDIM
           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,facezData,FACEZ)

           ! W velocity derivatives on external point:
           we   = 0.; dwdx = 0.; dwdy =0.;
           dwdz = 0.
           do i = 1 , ib_stencil

              ui = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              we   = we   + phile(i,1)*ui;
              dwdx = dwdx + phile(i,2)*ui;
              dwdy = dwdy + phile(i,3)*ui;
              dWdz = dwdz + phile(i,4)*ui;

           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

           end select

#endif
           
        end if

     enddo

  enddo

  ! 2nd point
  ! External Point Position:
  xbe2(IAXIS) = xp + nxp*enh2
  xbe2(JAXIS) = yp + nyp*enh2
#if NDIM == 3
  xbe2(KAXIS) = zp + nzp*enh2
#else
  xbe2(KAXIS) = 0.
#endif

  zpres2 = 0.
  zL2= 0.
  do presflag = CONSTANT_ZERO,CONSTANT_ONE

     ! N kij
#ifdef TANGENT_WITH_VORTICITY
        print*, 'TWO_NORMAL_POINTS does not work for TANGENT_WITH_VORTICITY CASE'
#else
     nkij = 1 + (1-presflag)*(NDIM-1) ! in 2D d/dx,d/dy; in 3D d/dx,d/dy,d/dz
#endif

     do gridind = 1,nkij

#ifdef TANGENT_WITH_VORTICITY
        print*, 'TWO_NORMAL_POINTS does not work for TANGENT_WITH_VORTICITY CASE'
        pause
#else
        ! Define Grids in case of vorticity and pressure:
        gridfl(:) = presflag*grdip(:) + (1-presflag)*grdu(:,gridind)
        ! Auxiliary deltas
        do idim = 1,NDIM
           delaux(idim) = (real(presflag)*dlip(idim)+real(1-presflag)*dlu(idim,gridind))*del(idim)
        enddo
#endif

        ! Obtain Stencil for External Point:
        call ib_stencils(xbe2,np,gridfl,del,coord,bsize,   & 
                         ielem(:,:,presflag+1),hl,COMPUTE_FORCES)

        ! Compute shape functions
        ! Positions of points on the stencil:
        xyz_stencil(1:ib_stencil,1:MDIM) = 0. 
        do idim = 1,NDIM
           xyz_stencil(1:ib_stencil,idim) = coord(idim) - 0.5*bsize(idim) + &
                real(ielem(1:ib_stencil,idim,presflag+1) - NGUARD - 1)*del(idim) + delaux(idim) 
        enddo

        ! Get interpolation functions:
        call ib_getInterpFunc(xbe2,xyz_stencil,del,derivflag,phile)

        if (presflag .eq. CONSTANT_ONE) then         ! Pressure

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,solnData,CENTER)

#ifndef TANGENT_WITH_VORTICITY
           dpdx2 = 0.; dpdy2 =0.;
           dpdz2 = 0.
           ppppxpx = 0.
           ppppxpy = 0.
           ppppypy = 0.

#endif
           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              p_i = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                      ielem(i,JAXIS,presflag+1), &
                                      ielem(i,KAXIS,presflag+1));  

              pixp = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                       ielem(i,JAXIS,presflag+1), &
                                       ielem(i,KAXIS,presflag+1));
              pixm = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                       ielem(i,JAXIS,presflag+1), &
                                       ielem(i,KAXIS,presflag+1));
              piyp = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                       ielem(i,JAXIS,presflag+1)+1, &
                                       ielem(i,KAXIS,presflag+1));
              piym = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                       ielem(i,JAXIS,presflag+1)-1, &
                                       ielem(i,KAXIS,presflag+1));

              pixpyp = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                         ielem(i,JAXIS,presflag+1)+1, &
                                         ielem(i,KAXIS,presflag+1));
              pixpym = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                         ielem(i,JAXIS,presflag+1)-1, &
                                         ielem(i,KAXIS,presflag+1));
              pixmyp = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                         ielem(i,JAXIS,presflag+1)+1, &
                                         ielem(i,KAXIS,presflag+1));
              pixmym = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                         ielem(i,JAXIS,presflag+1)-1, &
                                         ielem(i,KAXIS,presflag+1));

              zpres2 = zpres2 + phile(i,1)*p_i 
              dpdx2  = dpdx2  + phile(i,1)*(pixp-pixm)/(2.*dx)
              dpdy2  = dpdy2  + phile(i,1)*(piyp-piym)/(2.*dy)
              ppppxpx = ppppxpx + phile(i,1)*(pixp-2.*p_i+pixm)/(dx*dx)
              ppppypy = ppppypy + phile(i,1)*(piyp-2.*p_i+piym)/(dy*dy)
              ppppxpy = ppppxpy + phile(i,1)*(pixpyp-pixpym-pixmyp+pixmym)/(4.*dx*dy)
#ifndef TANGENT_WITH_VORTICITY
!              dpdx2 = dpdx2 + phile(i,2)*p_i;
!              dpdy2 = dpdy2 + phile(i,3)*p_i;
#if NDIM == MDIM
              dpdz2 = dpdz2 + phile(i,4)*p_i;
#endif             
#endif
           enddo

           ! Release Pointer
           call Grid_releaseBlkPtr(blockID,solnData,CENTER)

           ! Get Pressure approximation at surface marker: Acceleration +
           ! gravity effects.
           dpdn2 = -(     ubdd*nxp +      vbdd*nyp +      wbdd*nzp) + & ! -rho*Du/Dt * n 
                   (ins_gravX*nxp + ins_gravY*nyp + ins_gravZ*nzp);    ! +rho*    g * n
           zL2 = zpres2 - dpdn*enh2;

        else                                         ! Tangent stress

#ifdef TANGENT_WITH_VORTICITY
 
        print*, 'TWO_NORMAL_POINTS does not work for TANGENT_WITH_VORTICITY CASE'
#else

           select case(gridind)
           case(1)

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,facexData,FACEX)

           ! U velocity derivatives on external point:
           uen2 = 0.; dudx2 = 0.; dudy2 =0.;
           dudz2 = 0.
           ddudxdx2 = 0.
           ddudxdy2 = 0.
           ddudydy2 = 0.
           do i = 1 , ib_stencil

              ui = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              uixp = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              uixm = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              uiyp = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1)+1, &
                                           ielem(i,KAXIS,presflag+1));

              uiym = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1)-1, &
                                           ielem(i,KAXIS,presflag+1));

            uixpyp = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                           ielem(i,JAXIS,presflag+1)+1, &
                                           ielem(i,KAXIS,presflag+1));

            uixmyp = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                           ielem(i,JAXIS,presflag+1)+1, &
                                           ielem(i,KAXIS,presflag+1));

            uixpym = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                           ielem(i,JAXIS,presflag+1)-1, &
                                           ielem(i,KAXIS,presflag+1));

            uixmym = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                           ielem(i,JAXIS,presflag+1)-1, &
                                           ielem(i,KAXIS,presflag+1));

              uen2   = uen2   + phile(i,1)*ui
              dudx2  = dudx2  + phile(i,1)*(uixp-uixm)/(2.*dx)
              dudy2  = dudy2  + phile(i,1)*(uiyp-uiym)/(2.*dy)
              ddudxdx2 = ddudxdx2 + phile(i,1)*(uixp+uixm-2.*ui)/(dx*dx)
              ddudydy2 = ddudydy2 + phile(i,1)*(uiyp+uiym-2.*ui)/(dy*dy)
              ddudxdy2 = ddudxdy2 + & 
                      & phile(i,1)*(uixpyp-uixpym-uixmyp+uixmym)/(4.*dx*dy)

!              uen2   = uen2   + phile(i,1)*ui
!              dudx2  = dudx2  + phile(i,1)*(uixp-uixm)/(2.*dx)
!              dudy2  = dudy2  + phile(i,1)*(uiyp-uiym)/(2.*dy)

!              ue2   = ue2   + phile(i,1)*ui;
!              dudx2 = dudx2 + phile(i,2)*ui;
!              dudy2 = dudy2 + phile(i,3)*ui;
#if NDIM == MDIM
              dudz2 = dudz2 + phile(i,4)*ui;
#endif             
           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,facexData,FACEX)

           case(2)

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,faceyData,FACEY)

           ! V velocity derivatives on external point:
           ven2 = 0.; dvdx2 = 0.; dvdy2 =0.;
           dvdz2 = 0.
           ddvdxdx2 = 0.
           ddvdxdy2 = 0.
           ddvdydy2 = 0.
           do i = 1 , ib_stencil

              ui = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              uixp = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              uixm = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              uiyp = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1)+1, &
                                           ielem(i,KAXIS,presflag+1));

              uiym = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1)-1, &
                                           ielem(i,KAXIS,presflag+1));

            uixpyp = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                           ielem(i,JAXIS,presflag+1)+1, &
                                           ielem(i,KAXIS,presflag+1));

            uixmyp = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                           ielem(i,JAXIS,presflag+1)+1, &
                                           ielem(i,KAXIS,presflag+1));

            uixpym = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                           ielem(i,JAXIS,presflag+1)-1, &
                                           ielem(i,KAXIS,presflag+1));

            uixmym = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                           ielem(i,JAXIS,presflag+1)-1, &
                                           ielem(i,KAXIS,presflag+1));

              ven2   = ven2   + phile(i,1)*ui
              dvdx2  = dvdx2  + phile(i,1)*(uixp-uixm)/(2.*dx)
              dvdy2  = dvdy2  + phile(i,1)*(uiyp-uiym)/(2.*dy)
              ddvdxdx2 = ddvdxdx2 + phile(i,1)*(uixp+uixm-2.*ui)/(dx*dx)
              ddvdydy2 = ddvdydy2 + phile(i,1)*(uiyp+uiym-2.*ui)/(dy*dy)
              ddvdxdy2 = ddvdxdy2 + & 
                      & phile(i,1)*(uixpyp-uixpym-uixmyp+uixmym)/(4.*dx*dy)

!              ven2   = ven2   + phile(i,1)*ui
!              dvdx2  = dvdx2  + phile(i,1)*(uixp-uixm)/(2.*dx)
!              dvdy2  = dvdy2  + phile(i,1)*(uiyp-uiym)/(2.*dy)
!              ve2   = ve2   + phile(i,1)*ui; 
!              dvdx2 = dvdx2 + phile(i,2)*ui;
!              dvdy2 = dvdy2 + phile(i,3)*ui;
#if NDIM == MDIM
              dvdz2 = dvdz2 + phile(i,4)*ui;
#endif             
           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,faceyData,FACEY)


           case(3) 

#if NDIM == MDIM
           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,facezData,FACEZ)

           ! W velocity derivatives on external point:
           we2   = 0.; dwdx2 = 0.; dwdy2 =0.;
           dwdz2 = 0.
           do i = 1 , ib_stencil

              ui = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              we2   = we2   + phile(i,1)*ui;
              dwdx2 = dwdx2 + phile(i,2)*ui;
              dwdy2 = dwdy2 + phile(i,3)*ui;
              dWdz2 = dwdz2 + phile(i,4)*ui;

           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

           end select

#endif
           
        end if

     enddo

  enddo

!########################################################
  ! 3rd point
  ! External Point Position:
  xbe3(IAXIS) = xp + nxp*enh3
  xbe3(JAXIS) = yp + nyp*enh3
#if NDIM == 3
  xbe3(KAXIS) = zp + nzp*enh3
#else
  xbe3(KAXIS) = 0.
#endif

  do presflag = CONSTANT_ZERO,CONSTANT_ONE

     ! N kij
#ifdef TANGENT_WITH_VORTICITY
        print*, 'TWO_NORMAL_POINTS does not work for TANGENT_WITH_VORTICITY CASE'
#else
     nkij = 1 + (1-presflag)*(NDIM-1) ! in 2D d/dx,d/dy; in 3D d/dx,d/dy,d/dz
#endif

     do gridind = 1,nkij

#ifdef TANGENT_WITH_VORTICITY
        print*, 'TWO_NORMAL_POINTS does not work for TANGENT_WITH_VORTICITY CASE'
        pause
#else
        ! Define Grids in case of vorticity and pressure:
        gridfl(:) = presflag*grdip(:) + (1-presflag)*grdu(:,gridind)
        ! Auxiliary deltas
        do idim = 1,NDIM
           delaux(idim) = (real(presflag)*dlip(idim)+real(1-presflag)*dlu(idim,gridind))*del(idim)
        enddo
#endif

        ! Obtain Stencil for External Point:
        call ib_stencils(xbe3,np,gridfl,del,coord,bsize,   & 
                         ielem(:,:,presflag+1),hl,COMPUTE_FORCES)

        ! Compute shape functions
        ! Positions of points on the stencil:
        xyz_stencil(1:ib_stencil,1:MDIM) = 0. 
        do idim = 1,NDIM
           xyz_stencil(1:ib_stencil,idim) = coord(idim) - 0.5*bsize(idim) + &
                real(ielem(1:ib_stencil,idim,presflag+1) - NGUARD - 1)*del(idim) + delaux(idim) 
        enddo

        ! Get interpolation functions:
        call ib_getInterpFunc(xbe3,xyz_stencil,del,derivflag,phile)

        if (presflag .eq. CONSTANT_ONE) then         ! Pressure

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,solnData,CENTER)

#ifndef TANGENT_WITH_VORTICITY
!           dpdx2 = 0.; dpdy2 =0.;
!           dpdz2 = 0.
           ppppxpx3 = 0.
           ppppxpy3 = 0.
           ppppypy3 = 0.

#endif
           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              p_i = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                      ielem(i,JAXIS,presflag+1), &
                                      ielem(i,KAXIS,presflag+1));  

              pixp = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                       ielem(i,JAXIS,presflag+1), &
                                       ielem(i,KAXIS,presflag+1));
              pixm = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                       ielem(i,JAXIS,presflag+1), &
                                       ielem(i,KAXIS,presflag+1));
              piyp = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                       ielem(i,JAXIS,presflag+1)+1, &
                                       ielem(i,KAXIS,presflag+1));
              piym = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                       ielem(i,JAXIS,presflag+1)-1, &
                                       ielem(i,KAXIS,presflag+1));

              pixpyp = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                         ielem(i,JAXIS,presflag+1)+1, &
                                         ielem(i,KAXIS,presflag+1));
              pixpym = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                         ielem(i,JAXIS,presflag+1)-1, &
                                         ielem(i,KAXIS,presflag+1));
              pixmyp = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                         ielem(i,JAXIS,presflag+1)+1, &
                                         ielem(i,KAXIS,presflag+1));
              pixmym = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                         ielem(i,JAXIS,presflag+1)-1, &
                                         ielem(i,KAXIS,presflag+1));

              !zpres2 = zpres2 + phile(i,1)*p_i 
              !dpdx2  = dpdx2  + phile(i,1)*(pixp-pixm)/(2.*dx)
              !dpdy2  = dpdy2  + phile(i,1)*(piyp-piym)/(2.*dy)
              ppppxpx3 = ppppxpx3 + phile(i,1)*(pixp-2.*p_i+pixm)/(dx*dx)
              ppppypy3 = ppppypy3 + phile(i,1)*(piyp-2.*p_i+piym)/(dy*dy)
              ppppxpy3 = ppppxpy3 + phile(i,1)*(pixpyp-pixpym-pixmyp+pixmym)/(4.*dx*dy)
#ifndef TANGENT_WITH_VORTICITY
!              dpdx2 = dpdx2 + phile(i,2)*p_i;
!              dpdy2 = dpdy2 + phile(i,3)*p_i;
#if NDIM == MDIM
!              dpdz2 = dpdz2 + phile(i,4)*p_i;
#endif             
#endif
           enddo

           ! Release Pointer
           call Grid_releaseBlkPtr(blockID,solnData,CENTER)

           
        end if

     enddo

  enddo

!########################################################
  ! Project external and marker Velocity on the plane:
  ven = nxp*uen1 + nyp*ven1

  ! ves = ve - (ve*n) n
  ues  = uen1 - ven*nxp
  ves  = ven1 - ven*nyp
  wes  = 0. !wen1 - ven*nzp

  ! Project external and marker Velocity on the plane:
  ven = nxp*uen2 + nyp*ven2 !+ nzp*wen2

  ! ves = ve - (ve*n) n
  ues2  = uen2 - ven*nxp
  ves2  = ven2 - ven*nyp
  wes2  = 0. !wen2 - ven*nzp

  ! Linear Part:
  normt = 1./(enh2-enh1)
  dun_e = (ues2-ues)*normt
  dvn_e = (ves2-ves)*normt
  dwn_e = (wes2-wes)*normt

  ! vps = vp - (vp*n) n
  vpn  = nxp*ubd + nyp*vbd + nzp*wbd
  ups  = ubd - vpn*nxp
  vps  = vbd - vpn*nyp
  wps  = wbd - vpn*nzp

  ! Linear Part:
  normt = 1./enh1
  dun = (ues-ups)*normt
  dvn = (ves-vps)*normt
  dwn = (wes-wps)*normt

  ! Correction from diffusion equation (A. Posa):
  ! First versor in the local tangent velocity direction:
  if( (abs(dun)+abs(dvn)+abs(dwn)) .lt. 1.0e-14) then ! case zero velocity difference
  tx = 0.
  ty = 0.
  tz = 0.
  dpdxt = 0.
  else
  normt = 1./sqrt(dun**2. + dvn**2. + dwn**2.)
  tx = dun*normt
  ty = dvn*normt
  tz = dwn*normt
  endif

  ! slope
  if (NDIM==MDIM) pause 'slope not defined for 3D flows'
  sxxen1 = dudx
  sxyen1 = 0.5*(dudy + dvdx)
  syyen1 = dvdy
  syxen1 = 0.5*(dvdx + dudy)

  sxxen2 = dudx2
  sxyen2 = 0.5*(dudy2 + dvdx2)
  syyen2 = dvdy2
  syxen2 = 0.5*(dvdx2 + dudy2)

  ! Fvisc = Tau * n
  fvxen1 = nu*2.*(sxxen1*nxp+sxyen1*nyp)
  fvyen1 = nu*2.*(syxen1*nxp+syyen1*nyp)

  fvxen2 = nu*2.*(sxxen2*nxp+sxyen2*nyp)
  fvyen2 = nu*2.*(syxen2*nxp+syyen2*nyp)

  fvxw = fvxen1 - enh1/(enh2-enh1)*(fvxen2-fvxen1)
  fvyw = fvyen1 - enh1/(enh2-enh1)*(fvyen2-fvyen1)
  
  wz_2n_e1 = -dudy  + dvdx
  wz_2n_e2 = -dudy2 + dvdx2
  wz_2n_e  = wz_2n_e1 - enh1/(enh2-enh1)*(wz_2n_e2-wz_2n_e1) 

  wz_2n_lw = dvn*nxp - dun*nyp
  wz_2n_le = dvn_e*nxp - dun_e*nyp

  fvxen1_c = -nu*enh1*(nxp*nxp*ddudxdx+nxp*nyp*ddudxdy*2.+nyp*nyp*ddudydy)
  fvyen1_c = -nu*enh1*(nxp*nxp*ddvdxdx+nxp*nyp*ddvdxdy*2.+nyp*nyp*ddvdydy)

  fvxw_c = fvxen1 + fvxen1_c
  fvyw_c = fvyen1 + fvyen1_c

  wz_2n_e1_c = fvyen1_c/nu*nxp - fvxen1_c/nu*nyp
  wz_2n_w_c  = wz_2n_e1 + wz_2n_e1_c
 
  fvxen2_c = -nu*enh2*(nxp*nxp*ddudxdx2+nxp*nyp*ddudxdy2*2.+nyp*nyp*ddudydy2)
  fvyen2_c = -nu*enh2*(nxp*nxp*ddvdxdx2+nxp*nyp*ddvdxdy2*2.+nyp*nyp*ddvdydy2)

  fvxen2_c2 = ((fvxen2_c*tx+fvyen2_c*ty) - (fvxen1_c*tx+fvyen1_c*ty))/(enh2-enh1)*tx*enh2*enh2*0.5
  fvyen2_c2 = ((fvxen2_c*tx+fvyen2_c*ty) - (fvxen1_c*tx+fvyen1_c*ty))/(enh2-enh1)*ty*enh2*enh2*0.5
               
  fvxw_c2 = fvxen2 + fvxen2_c
  fvyw_c2 = fvyen2 + fvyen2_c

  wz_2n_e2_c = fvyen2_c/nu*nxp - fvxen2_c/nu*nyp
  wz_2n_w_c2  = wz_2n_e2 + wz_2n_e2_c
 
  wz_2n_e2_c2 = fvyen2_c2/nu*nxp - fvxen2_c2/nu*nyp

  pppxn1 = dpdx*nxp + dpdy*nyp
  pppxn2 = dpdx2*nxp + dpdy2*nyp
  pres_c2 = zpres2 - pppxn2*enh2

  pppxt1 = dpdx*tx + dpdy*ty
  pppxt2 = dpdx2*tx + dpdy2*ty
  ppppxnpxt2 = nxp*ppppxpx*tx + nxp*ppppxpy*ty + nyp*ppppxpy*tx + nyp*ppppypy*ty

  ppppxnpxt3 = (nxp*ppppxpx3*tx + nxp*ppppxpy3*ty + nyp*ppppxpy3*tx + nyp*ppppypy3*ty)

  fvxen2_pc = -enh2*pppxt2*tx
  fvyen2_pc = -enh2*pppxt2*ty

  wz_2n_e2_pc = fvyen2_pc/nu*nxp - fvxen2_pc/nu*nyp
  
  fvxen2_ppc = -enh2*enh2*ppppxnpxt2*tx*0.5
  fvyen2_ppc = -enh2*enh2*ppppxnpxt2*ty*0.5

  wz_2n_e2_ppc = fvyen2_ppc/nu*nxp - fvxen2_ppc/nu*nyp
  
  fvxen2_ppc3 = -enh3*enh3*ppppxnpxt3*tx*0.5/enh3*enh2
  fvyen2_ppc3 = -enh3*enh3*ppppxnpxt3*ty*0.5/enh3*enh2

  wz_2n_e2_ppc3 = fvyen2_ppc3/nu*nxp - fvxen2_ppc3/nu*nyp

!#ifdef TEST_COMPARE
  write(9200+gr_meshMe, '(30f20.12)') particleData(GLOB_PART_PROP), pres_c2, zL2, zpres2
  write(9300+gr_meshMe, '(30f20.12)') particleData(GLOB_PART_PROP), wz_2n_e2+wz_2n_e2_c-wz_2n_e2_ppc, &
                             wz_2n_e2, wz_2n_e2_c, -wz_2n_e2_ppc

  write(9500+gr_meshMe, '(30f20.8)') particleData(GLOB_PART_PROP), zL, zL2,  &
                       zpres, zpres2, dpdx, dpdy, dpdx2, dpdy2,              &
                       ppppxpx, ppppxpy, ppppypy, pppxn2, pppxn1, pppxt2, pppxt1, &
                       ppppxnpxt2, ppppxnpxt3, nxp, nyp, tx, ty 

  write(9600+gr_meshMe, '(30f20.8)') particleData(GLOB_PART_PROP), zL, zL2,  &
                       zpres, zpres2,                                        &
                       fvxw, fvyw, fvxen1, fvyen1, fvxen2, fvyen2,           &
                       wz_2n_e, wz_2n_e1, wz_2n_e2, wz_2n_lw, wz_2n_le

  write(9700+gr_meshMe, '(30f20.8)') particleData(GLOB_PART_PROP), zL, zL2,  &
                       zpres, zpres2,                                        &
                       fvxw_c, fvyw_c, fvxen1, fvyen1, fvxen1_c, fvyen1_c,           &
                       wz_2n_w_c, wz_2n_e1, wz_2n_e1_c, wz_2n_lw, wz_2n_le

  write(9800+gr_meshMe, '(30f20.8)') particleData(GLOB_PART_PROP), zL, zL2,  &
                       zpres, zpres2,                                        &
                       fvxw_c2, fvyw_c2, fvxen2, fvyen2, fvxen2_c, fvyen2_c,           &
                       wz_2n_w_c2, wz_2n_e2, wz_2n_e2_c, wz_2n_lw, wz_2n_le, &
                       pres_c2, -pppxn2*enh2
  
  write(9900+gr_meshMe, '(30f20.8)') particleData(GLOB_PART_PROP),           &
                       wz_2n_w_c2, wz_2n_e2, wz_2n_e2_c,                     &
                       wz_2n_w_c,  wz_2n_e1, wz_2n_e1_c,                     &
                       wz_2n_e2_pc, wz_2n_e2_ppc, wz_2n_e2_c2, wz_2n_e2_ppc3
!#endif

#endif /* TWO_NORMAL_POINTS_EUL */

!#######################################################################################################################
!#######################################################################################################################
  ! Assign pressure and viscous forces to particleData:
#ifdef INS_CONSTDENS
  ! Case constant density:
  particleData(PRES_PART_PROP) = zL !+ ins_gravX*(xp-gr_imin) + &
                                    !  ins_gravY*(yp-gr_jmin) + &
                                    !  ins_gravZ*(zp-gr_kmin)
  particleData(PEX0_PART_PROP) = particleData(PEXT_PART_PROP)
  particleData(PEXT_PART_PROP) = zpres !+ ins_gravX*(xbe(IAXIS)-gr_imin) + &
                                       !  ins_gravY*(xbe(JAXIS)-gr_jmin) + &
                                       !  ins_gravZ*(xbe(KAXIS)-gr_kmin)


#ifdef TANGENT_WITH_VORTICITY 

  ! Fvisc = -nu (N x W)
  particleData(FXVI_PART_PROP) = (nzp*nuwy-nyp*nuwz)
  particleData(FYVI_PART_PROP) = (nxp*nuwz-nzp*nuwx)
  particleData(FZVI_PART_PROP) = (nyp*nuwx-nxp*nuwy)

#else

#ifdef USE_CF

  ! Project external and marker Velocity on the plane:
  ven = nxp*ue + nyp*ve + nzp*we

  ! ves = ve - (ve*n) n
  ues  = ue - ven*nxp
  ves  = ve - ven*nyp
  wes  = we - ven*nzp

  ! vps = vp - (vp*n) n
  vpn  = nxp*ubd + nyp*vbd + nzp*wbd
  ups  = ubd - vpn*nxp
  vps  = vbd - vpn*nyp
  wps  = wbd - vpn*nzp

  ! Linear Part:
  normt = 1./h
  dun = (ues-ups)*normt
  dvn = (ves-vps)*normt
  dwn = (wes-wps)*normt

  ! Correction from diffusion equation (A. Posa):
  ! First versor in the local tangent velocity direction:
  if( (abs(dun)+abs(dvn)+abs(dwn)) .lt. 1.0e-14) then ! case zero velocity difference
  tx = 0.
  ty = 0.
  tz = 0.
  dpdxt = 0.
  else
  normt = 1./sqrt(dun**2. + dvn**2. + dwn**2.)
  tx = dun*normt
  ty = dvn*normt
  tz = dwn*normt

  ! Compute dpdxt from pressure gradients directly, use only hydrodynamic dpdxt:
  dpdxt = (dpdx-ins_gravX)*tx + (dpdy-ins_gravY)*ty + (dpdz-ins_gravZ)*tz


!  x = xbe(1)
!  y = xbe(2)
!  write(10000+gr_meshMe, '(20f20.8)') particleData(GLOB_PART_PROP), xbe(1), xbe(2), &
!                 nxp, nyp, tx, ty, zpres, ue, ve, dpdx, dpdy, dudx, dudy, dvdx, dvdy, dpdxt

!  write(11000+gr_meshMe, '(20f20.8)') particleData(GLOB_PART_PROP), xbe(1), xbe(2), &
!                 nxp, nyp, tx, ty, -0.25*(cos(2.*pi*x) + cos(2.*pi*y) ), &
!                 -cos(pi*x)*sin(pi*y), sin(pi*x)*cos(pi*y),              &
!                 0.25*2.*pi*sin(2.*pi*x), 0.25*2.*pi*sin(2.*pi*y),       &
!                 pi*sin(pi*x)*sin(pi*y), -pi*cos(pi*x)*cos(pi*y),        &
!                 pi*cos(pi*x)*cos(pi*y), -pi*sin(pi*x)*sin(pi*y),        &
!                 pi*sin(pi*x)*sin(pi*y)*tx-pi*sin(pi*x)*sin(pi*y)*ty 

   dun_l_2p = dun
   dvn_l_2p = dvn

#ifdef TWO_POINTSP
  ! Compute dpdxt from other pressure values:
  ! External Point Position:
  normt = h/2.
  xbe2(IAXIS) = xp + nxp*h + tx*normt
  xbe2(JAXIS) = yp + nyp*h + ty*normt
#if NDIM == 3
  xbe2(KAXIS) = zp + nzp*h + tz*normt
#else
  xbe2(KAXIS) = 0.
#endif
  zpres2 = 0.
  presflag = CONSTANT_ONE
  gridind  = 1
  ! Define Grids in case of pressure:
  gridfl(:) = grdip(:)
  ! Auxiliary deltas
  do idim = 1,MDIM
     delaux(idim) = (dlip(idim))*del(idim)
  enddo

  !! Point 2:
  ! Obtain Stencil for External Point:
  call ib_stencils(xbe2,np,gridfl,del,coord,bsize,   &
                   ielem(:,:,presflag+1),hl,COMPUTE_FORCES)
  ! Compute shape functions
  ! Positions of points on the stencil:
  xyz_stencil(1:ib_stencil,1:MDIM) = 0.
  do idim = 1,NDIM
     xyz_stencil(1:ib_stencil,idim) = coord(idim) - 0.5*bsize(idim) + &
                real(ielem(1:ib_stencil,idim,presflag+1) - NGUARD - 1)*del(idim) + delaux(idim)
  enddo
  ! Get interpolation functions:
  call ib_getInterpFunc(xbe2,xyz_stencil,del,0,phile)
  ! Point to cell centered Variables:
  call Grid_getBlkPtr(blockID,solnData,CENTER)
  ! Value of the function in xbe:
  do i = 1 , ib_stencil
      p_i = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                              ielem(i,JAXIS,presflag+1), &
                              ielem(i,KAXIS,presflag+1));
      zpres2 = zpres2 + phile(i,1)*p_i
  enddo
  ! Release Pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  ! Substract hydrostatic pressure to total pressure in point 2:
  zpres2 = zpres2 - ins_gravX*(xbe2(IAXIS)-gr_imin) - &
                    ins_gravY*(xbe2(JAXIS)-gr_jmin) - &
                    ins_gravZ*(xbe2(KAXIS)-gr_kmin)

  !! Point 3:
  xbe3(IAXIS) = xp + nxp*h - tx*normt
  xbe3(JAXIS) = yp + nyp*h - ty*normt
#if NDIM == 3
  xbe3(KAXIS) = zp + nzp*h - tz*normt
#else
  xbe3(KAXIS) = 0.
#endif
  zpres3 = 0.
  ! Obtain Stencil for External Point:
  call ib_stencils(xbe3,np,gridfl,del,coord,bsize,   &
                   ielem(:,:,presflag+1),hl,COMPUTE_FORCES)
  ! Compute shape functions
  ! Positions of points on the stencil:
  xyz_stencil(1:ib_stencil,1:MDIM) = 0.
  do idim = 1,NDIM
     xyz_stencil(1:ib_stencil,idim) = coord(idim) - 0.5*bsize(idim) + &
                real(ielem(1:ib_stencil,idim,presflag+1) - NGUARD - 1)*del(idim) + delaux(idim)
  enddo
  ! Get interpolation functions:
  call ib_getInterpFunc(xbe3,xyz_stencil,del,0,phile)
  ! Point to cell centered Variables:
  call Grid_getBlkPtr(blockID,solnData,CENTER)
  ! Value of the function in xbe:
  do i = 1 , ib_stencil
      p_i = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                              ielem(i,JAXIS,presflag+1), &
                              ielem(i,KAXIS,presflag+1));
      zpres3 = zpres3 + phile(i,1)*p_i
  enddo
  ! Release Pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  ! Substract hydrostatic pressure to total pressure in point 3:
  zpres3 = zpres3 - ins_gravX*(xbe3(IAXIS)-gr_imin) - &
                    ins_gravY*(xbe3(JAXIS)-gr_jmin) - &
                    ins_gravZ*(xbe3(KAXIS)-gr_kmin)

  ! Finally dpdxt:
  !dpdxt = (zpres2-zpres)/normt ! This 1st order approx would save us from
                                ! computing zpres3
  dpdxt = (zpres2-zpres3)/(2.*normt) 

#endif /* TWO POINT flag */

#ifdef PRES_GRAD_CORR

  If(NDIM == MDIM) then
    print*, 'The PRES_GRAD_CORR does not work for 3D flows now'
    stop
  Endif

  pppx = 0.
  pppy = 0.
  pppxpx = 0.
  pppxpy = 0.
  pppypx = 0.
  pppypy = 0.

  presflag = CONSTANT_ONE
  gridind  = 1
  ! Define Grids in case of pressure:
  gridfl(:) = grdip(:)
  ! Auxiliary deltas
  do idim = 1,MDIM
     delaux(idim) = (dlip(idim))*del(idim)
  enddo

  ! Obtain Stencil for External Point:
  call ib_stencils(xbe,np,gridfl,del,coord,bsize,   &
                   ielem(:,:,presflag+1),hl,COMPUTE_FORCES)

  ! Compute shape functions
  ! Positions of points on the stencil:
  xyz_stencil(1:ib_stencil,1:MDIM) = 0.
  do idim = 1,NDIM
     xyz_stencil(1:ib_stencil,idim) = coord(idim) - 0.5*bsize(idim) + &
                real(ielem(1:ib_stencil,idim,presflag+1) - NGUARD - 1)*del(idim) + delaux(idim)
  enddo
  ! Get interpolation functions:
  call ib_getInterpFunc(xbe,xyz_stencil,del,derivflag,phile)
  ! Point to cell centered Variables:
  call Grid_getBlkPtr(blockID,solnData,CENTER)
  ! Value of the function in xbe:
  do i = 1 , ib_stencil
      p_i = solnData(PPPX_VAR,ielem(i,IAXIS,presflag+1), &
                              ielem(i,JAXIS,presflag+1), &
                              ielem(i,KAXIS,presflag+1));
      pppx    = pppx   + phile(i,1)*p_i
      pppxpx  = pppxpx + phile(i,2)*p_i
      pppxpy  = pppxpy + phile(i,3)*p_i

      p_i = solnData(PPPY_VAR,ielem(i,IAXIS,presflag+1), &
                              ielem(i,JAXIS,presflag+1), &
                              ielem(i,KAXIS,presflag+1));
      pppy    = pppy   + phile(i,1)*p_i
      pppypx  = pppypx + phile(i,2)*p_i
      pppypy  = pppypy + phile(i,3)*p_i
  enddo

!  print*, 'ib_stencil', ib_stencil
!  print*, 'stencil', phile
!  print*, 'xbe', xbe, xp, yp
!  print*, 'pppx', pppx, xbe(1)+xbe(2)**2, pppxpx, 1.0, pppxpy, xbe(2)
!  print*, 'pppy', pppy, sin(xbe(1))+cos(xbe(2)), pppypx, cos(xbe(1)),  pppypy, -sin(xbe(2))
 
!  write(9000+gr_meshMe, '(10f20.8)'), pppx, xbe(1)+xbe(2)**2, pppxpx, 1.0, pppxpy, 2*xbe(2)
!  write(8000+gr_meshMe, '(10f20.8)'), pppy, sin(xbe(1))+cos(xbe(2)), pppypx, cos(xbe(1)),  pppypy, -sin(xbe(2))

  ! Release Pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  ! Get the dpdxt 
  dpdxt = (pppx-ins_gravX)*tx + (pppy-ins_gravY)*ty

!  write(6000+gr_meshMe, '(10f20.8)') particleData(GLOB_PART_PROP), xp, yp, dpdxt

!  dpdxtdxn = nxp*pppxpx*tx + nxp*pppxpy*ty + &
!           & nyp*pppypx*tx + nyp*pppypy*ty + &
!           & nxp*tx*pppx   + nyp*ty*pppy 
  dpdxtdxn = nxp*pppxpx*tx + nxp*pppxpy*ty + &
           & nyp*pppypx*tx + nyp*pppypy*ty

  dpdxt  = dpdxt - dpdxtdxn*h

!  x = xbe(1)
!  y = xbe(2)
!  write(30000+gr_meshMe, '(20f20.8)') particleData(GLOB_PART_PROP), xbe(1), xbe(2), &
!                 nxp, nyp, tx, ty, pppx, pppy, pppxpx, pppxpy, pppypx, pppypy
  
  ! p = -0.25*(cos(2.*pi*x) + cos(2.*pi*y) )
  ! u = -cos(pi*x)*sin(pi*y),
  ! v = sin(pi*x)*cos(pi*y)
!  write(31000+gr_meshMe, '(20f20.8)') particleData(GLOB_PART_PROP), xbe(1), xbe(2), &
!                 nxp, nyp, tx, ty,                                                  &
!                 0.25*2.*pi*sin(2.*pi*x),       &
!                 0.25*2.*pi*sin(2.*pi*y),       &
!                 0.25*2.*pi*2.*pi*cos(2.*pi*x), &
!                 0.0,                           &
!                 0.0,                           &
!                 0.25*2.*pi*2.*pi*cos(2.*pi*y)

#endif  /* PRES_GRAD_CORR */


#ifdef PRES_BDL_EQ

  If(NDIM == MDIM) then
    print*, 'The PRES_BDL_EQ does not work for 3D flows now'
    stop
  Endif

  pupx = 0.
  pupy = 0.
  pupxpx = 0.
  pupxpy = 0.
  pupypx = 0.
  pupypy = 0.

  pvpx = 0.
  pvpy = 0.
  pvpxpx = 0.
  pvpxpy = 0.
  pvpypx = 0.
  pvpypy = 0.

  presflag = CONSTANT_ZERO
  nkij  = 1 + (1-presflag)*(NDIM-1)


  Do gridind = 1, nkij

  ! Define Grids in case of pressure:
  gridfl(:) = grdu(:,gridind)
  ! Auxiliary deltas
  do idim = 1,MDIM
     delaux(idim) = dlu(idim,gridind)*del(idim)
  enddo

  ! Obtain Stencil for External Point:
  call ib_stencils(xbe,np,gridfl,del,coord,bsize,   &
                   ielem(:,:,presflag+1),hl,COMPUTE_FORCES)

  ! Compute shape functions
  ! Positions of points on the stencil:
  xyz_stencil(1:ib_stencil,1:MDIM) = 0.
  do idim = 1,NDIM
     xyz_stencil(1:ib_stencil,idim) = coord(idim) - 0.5*bsize(idim) + &
                real(ielem(1:ib_stencil,idim,presflag+1) - NGUARD - 1)*del(idim) + delaux(idim)
  enddo
  ! Get interpolation functions:
  call ib_getInterpFunc(xbe,xyz_stencil,del,derivflag,phile)
  ! Point to cell centered Variables:

           select case(gridind)
           case(1)

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,facexData,FACEX)

           ! U velocity derivatives on external point:
           do i = 1 , ib_stencil

              ui = facexData(PPX_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              pupx   = pupx   + phile(i,1)*ui;
              pupxpx = pupxpx + phile(i,2)*ui;
              pupxpy = pupxpy + phile(i,3)*ui;

              ui = facexData(PPY_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              pupy   = pupy   + phile(i,1)*ui;
              pupypx = pupypx + phile(i,2)*ui;
              pupypy = pupypy + phile(i,3)*ui;
           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,facexData,FACEX)

           case(2)

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,faceyData,FACEY)

           ! V velocity derivatives on external point:
           do i = 1 , ib_stencil

              ui = faceyData(PPX_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              pvpx   = pvpx   + phile(i,1)*ui; 
              pvpxpx = pvpxpx + phile(i,2)*ui;
              pvpxpy = pvpxpy + phile(i,3)*ui;
              
              ui = faceyData(PPY_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              pvpy   = pvpy   + phile(i,1)*ui; 
              pvpypx = pvpypx + phile(i,2)*ui;
              pvpypy = pvpypy + phile(i,3)*ui;
           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

           end select

  Enddo

  ! Get the dpdxt from boundary layer equation along local streamline
  ! dpdxt = rhou*(-Du/Dt+nu*partial^2_ut/partial_xn^2 + gt) 
  dpdxt = -(ue*tx*pupx + ue*ty*pvpx + ve*tx*pupy + ve*ty*pvpy) +  &
          nu*tx*(nxp*nxp*pupxpx+2.*nxp*nyp*pupxpy+nyp*nyp*pupypy) +     &
          nu*ty*(nxp*nxp*pvpxpx+2.*nxp*nyp*pvpxpy+nyp*nyp*pvpypy) +     &
          (ins_gravX)*tx + (ins_gravY)*ty

#ifdef PRES_NS_EQ
  ! Get the dpdxt from NS equations along local streamline
  ! dpdxt = rhou*(-Du/Dt+nu*partial^2_ut/partial_xn^2 +
  ! nu*partial^2_ut/partial_xt^2 + gt) 
  dpdxt = -(ue*tx*pupx + ue*ty*pvpx + ve*tx*pupy + ve*ty*pvpy) +  &
          nu*tx*(nxp*nxp*pupxpx+2.*nxp*nyp*pupxpy+nyp*nyp*pupypy) +     &
          nu*ty*(nxp*nxp*pvpxpx+2.*nxp*nyp*pvpxpy+nyp*nyp*pvpypy) +     &
          nu*tx*(tx*tx*pupxpx+2.*tx*ty*pupxpy+ty*ty*pupypy) +     &
          nu*ty*(ty*ty*pvpxpx+2.*ty*ty*pvpxpy+ty*ty*pvpypy) +     &
          (ins_gravX)*tx + (ins_gravY)*ty

!  x = xbe(1)
!  y = xbe(2)
!  write(20000+gr_meshMe, '(20f20.8)') particleData(GLOB_PART_PROP), xbe(1), xbe(2), &
!                 nxp, nyp, tx, ty, pupx, pupy, pupxpx, pupxpy, pupypx, pupypy,      & 
!                                   pvpx, pvpy, pvpxpx, pvpxpy, pvpypx, pvpypy,      &
!                                   dpdxt
  ! p = -0.25*(cos(2.*pi*x) + cos(2.*pi*y) )
  ! u = -cos(pi*x)*sin(pi*y),
  ! v = sin(pi*x)*cos(pi*y)
!  write(21000+gr_meshMe, '(20f20.8)') particleData(GLOB_PART_PROP), xbe(1), xbe(2), &
!                 nxp, nyp, tx, ty,                                                  &
!                 pi*sin(pi*x)*sin(pi*y),        &
!                -pi*cos(pi*x)*cos(pi*y),        &
!                 pi*pi*cos(pi*x)*sin(pi*y),     &
!                 pi*pi*sin(pi*x)*cos(pi*y),     &
!                 pi*pi*sin(pi*x)*cos(pi*y),     &
!                 pi*pi*cos(pi*x)*sin(pi*y),     &                 
!                 pi*cos(pi*x)*cos(pi*y),        &
!                -pi*sin(pi*x)*sin(pi*y),        &
!                -pi*pi*sin(pi*x)*cos(pi*y),     &
!                -pi*pi*cos(pi*x)*sin(pi*y),     &
!                -pi*pi*cos(pi*x)*sin(pi*y),     &
!                -pi*pi*sin(pi*x)*cos(pi*y),     & 
!                 0.25*2.*pi*sin(2.*pi*x)*tx+0.25*2.*pi*cos(2.*pi*y)*ty 

#endif
#endif /* PRES_BDL_EQ  */

  endif ! case zero velocity difference

#ifdef TEST_COMPARE
  write(7000+gr_meshMe, '(10f20.8)') particleData(GLOB_PART_PROP), xp, yp, dpdxt, & 
                    zL, nxp, nyp, tx, ty, nxp*tx+nyp*ty 
#endif

  ! Diffusion correction - 1/rho*dp/dxt*h/(2nu) * t:
  normt = h/(2.*nu)
  !normt = h/(1.*nu) ! Shizhao for keeping the gradient at the external point

  ! Using Strain tensor on the external point:
  ! Strain velocities tensor:
  exx = dudx
  exy = 0.5*(dudy + dvdx)
  exz = 0.5*(dudz + dwdx)
  eyy = dvdy
  eyz = 0.5*(dvdz + dwdy)
  ezz = dwdz

  ! Fvisc = Tau * n
  dun = 2.*(exx*nxp+exy*nyp+exz*nzp) - (dpdxt*tx)*normt
  dvn = 2.*(exy*nxp+eyy*nyp+eyz*nzp) - (dpdxt*ty)*normt
  dwn = 2.*(exz*nxp+eyz*nyp+ezz*nzp) - (dpdxt*tz)*normt

  ! Shizhao test the pressure correction
  !dun = 2.*(exx*nxp+exy*nyp+exz*nzp)
  !dvn = 2.*(exy*nxp+eyy*nyp+eyz*nzp)
  !dwn = 2.*(exz*nxp+eyz*nyp+ezz*nzp)

  ! Using Linear du/dxn:
!  dun = dun - (dpdxt*tx)*normt
!  dvn = dvn - (dpdxt*ty)*normt
!  dwn = dwn - (dpdxt*tz)*normt
 
  ! With quadratic correction:
  !dune = dudx*nxp + dudy*nyp + dudz*nzp
  !dvne = dvdx*nxp + dvdy*nyp + dvdz*nzp
  !dwne = dwdx*nxp + dwdy*nyp + dwdz*nzp 
  !dun = 2.*(ues-ups)/h - dune
  !dvn = 2.*(ves-vps)/h - dvne
  !dwn = 2.*(wes-wps)/h - dwne

  ! Fvisc = nu dv/dn - 1/rho*dp/dxt*h/2 * t, here rho=1:
  particleData(FXVI_PART_PROP) = nu*dun 
  particleData(FYVI_PART_PROP) = nu*dvn 
  particleData(FZVI_PART_PROP) = nu*dwn 
 
  ! Vorticity:
  ! In Z dir: wz = dv/dx - du/dy
  zv(1) = dvn*nxp - dun*nyp
  ! In X dir: wx = dw/dy - dv/dz
  zv(2) = dwn*nyp - dvn*nzp
  ! In Y dir: wy = du/dz - dw/dx
  zv(3) = dun*nzp - dwn*nxp

#ifdef TWO_NORMAL_POINTS
  dun_l_vg = 2.*(exx*nxp+exy*nyp+exz*nzp)
  dvn_l_vg = 2.*(exy*nxp+eyy*nyp+eyz*nzp)
   
  dun_c = -(dpdxt*tx)*normt
  dvn_c = -(dpdxt*ty)*normt
  dun_c2 = -2.0*(dpdxt*tx)*normt
  dvn_c2 = -2.0*(dpdxt*ty)*normt

  dun_vgc  = dun_l_vg + dun_c
  dvn_vgc  = dvn_l_vg + dvn_c
  dun_vgc2 = dun_l_vg + dun_c2
  dvn_vgc2 = dvn_l_vg + dvn_c2

  dun_2pc  = dun_l_2p + dun_c
  dvn_2pc  = dvn_l_2p + dvn_c
  dun_2pc2 = dun_l_2p + dun_c2
  dvn_2pc2 = dvn_l_2p + dvn_c2

  wz_l_2p = dvn_l_2p*nxp - dun_l_2p*nyp
  wz_l_vg = dvn_l_vg*nxp - dun_l_vg*nyp
  wz_c    = dvn_c*nxp    - dun_c*nyp
  wz_c2   = dvn_c2*nxp   - dun_c2*nyp
  wz_vgc  = dvn_vgc*nxp  - dun_vgc*nyp
  wz_vgc2 = dvn_vgc2*nxp - dun_vgc2*nyp
  wz_2pc  = dvn_2pc*nxp  - dun_2pc*nyp
  wz_2pc2 = dvn_2pc2*nxp - dun_2pc2*nyp
#endif

#else /* USE_CF */


  ! Strain velocities tensor:
  exx = dudx
  exy = 0.5*(dudy + dvdx)
  exz = 0.5*(dudz + dwdx)
  eyy = dvdy 
  eyz = 0.5*(dvdz + dwdy)
  ezz = dwdz 

  ! Fvisc = Tau * n
  particleData(FXVI_PART_PROP) = 2.*nu*(exx*nxp+exy*nyp+exz*nzp)
  particleData(FYVI_PART_PROP) = 2.*nu*(exy*nxp+eyy*nyp+eyz*nzp)
  particleData(FZVI_PART_PROP) = 2.*nu*(exz*nxp+eyz*nyp+ezz*nzp)

  ! Vorticity:
  ! In Z dir: wz = dv/dx - du/dy
  zv(1) = dvdx - dudy
  ! In X dir: wx = dw/dy - dv/dz
  zv(2) = dwdy - dvdz
  ! In Y dir: wy = du/dz - dw/dx
  zv(3) = dudz - dwdx


#endif
#endif

#else /* not INS_CONSTDENS */

  call Driver_abortFlash("Variable density particle pressure and tangent stress not defined.")

#endif

  ! Vorticity components:
  particleData(VORZ_PART_PROP) = zv(1)
  particleData(VORX_PART_PROP) = zv(2)
  particleData(VORY_PART_PROP) = zv(3)

#ifdef TWO_NORMAL_POINTS
  particleData(FXVI_PART_PROP) = nu*dun_2n_vg 
  particleData(FYVI_PART_PROP) = nu*dvn_2n_vg 
#endif
#ifdef TWO_NORMAL_POINTS_EUL
  particleData(FXVI_PART_PROP) = fvxw_c2
  particleData(FYVI_PART_PROP) = fvyw_c2 
  particleData(PRES_PART_PROP) = pres_c2
#endif

#ifdef TEST_COMPARE
  write(9000+gr_meshMe, '(30f20.8)') particleData(GLOB_PART_PROP), xp, yp, & 
                    nxp, nyp, tx, ty, zL, nu*dun, nu*dvn, zv(1), dpdxt,    &
                    nu*2.*(exx*nxp+exy*nyp+exz*nzp), -nu*(dpdxt*tx)*normt,       &
                    nu*2.*(exy*nxp+eyy*nyp+eyz*nzp), -nu*(dpdxt*ty)*normt, & 
                    nu*(ues-ups)/h, nu*(ves-vps)/h

  write(9100+gr_meshMe, '(30f20.8)') particleData(GLOB_PART_PROP), zL, &
                       nu*dun_vgc, nu*dvn_vgc, wz_vgc,wz_l_vg, wz_c,  dpdxt,          &
                       nu*dun_l_vg, nu*dvn_l_vg,                       &
                       nu*dun_c,    nu*dvn_c

  write(9200+gr_meshMe, '(30f20.8)') particleData(GLOB_PART_PROP), zL, &
                       nu*dun_vgc2, nu*dvn_vgc2, wz_vgc2,wz_l_vg, wz_c2,  dpdxt,       &
                       nu*dun_l_vg, nu*dvn_l_vg,                       &
                       nu*dun_c2,   nu*dvn_c2

  write(9300+gr_meshMe, '(30f20.8)') particleData(GLOB_PART_PROP), zL, &
                       nu*dun_2pc, nu*dvn_2pc, wz_2pc,wz_l_2p, wz_c,  dpdxt,          &
                       nu*dun_l_2p, nu*dvn_l_2p,                       &
                       nu*dun_c,    nu*dvn_c

  write(9400+gr_meshMe, '(30f20.8)') particleData(GLOB_PART_PROP), zL, &
                       nu*dun_2pc2, nu*dvn_2pc2, wz_2pc2, wz_l_2p, wz_c2,  dpdxt,       &
                       nu*dun_l_2p, nu*dvn_l_2p,                       &
                       nu*dun_c2,   nu*dvn_c2

#endif

  return

end subroutine ib_distributedForces


