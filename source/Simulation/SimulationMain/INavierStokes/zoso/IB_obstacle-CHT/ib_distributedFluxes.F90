

#include "constants.h"
#include "Flash.h"
#include "ImBound.h"

subroutine ib_distributedFluxes(blockID, particleData, vortx, vorty, vortz)

  use Grid_Data, only : gr_meshMe

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,      &
                             Grid_getDeltas, Grid_getBlkCenterCoords, &
                             Grid_getBlkPhysicalSize

  use ImBound_Data, only : ib_nu,ib_stencil,ib_alphax,ib_alphay,ib_alphaz,ib_dt

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

  integer :: presflag,gridfl(MDIM),i,idim,gridind,nkij
  integer :: ielem(ib_stencil,MDIM,CONSTANT_TWO)
  real :: dpdn,eps
  real :: delaux(MDIM),xyz_stencil(ib_stencil,MDIM),phile(ib_stencil,NDIM+1)

  real :: zpres,p_i,zv(MDIM),nuwx,nuwy,nuwz

  real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData,facezData

  integer, parameter :: derivflag = 1 ! Give interpolation functions and their derivatives

  real ::nu, dt

  nu = ib_nu
  dt = ib_dt

  ue =0.; ve=0.; we=0.; 

  dudx=0.; dudy=0.; dudz=0.;
  dvdx=0.; dvdy=0.; dvdz=0.; 
  dwdx=0.; dwdy=0.; dwdz=0.;

  exx=0;  exy=0.; exz=0.;
  eyy=0.; eyz=0.;
  ezz=0.;

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
  xbe(IAXIS) = xp
  xbe(JAXIS) = yp
#if NDIM == 3
  xbe(KAXIS) = zp
#else
  xbe(KAXIS) = 0.
#endif

  xbe2(IAXIS) = xp + nxp*1.5*dx
  xbe2(JAXIS) = yp + nyp*1.5*dx
#if NDIM == 3
  xbe2(KAXIS) = zp + nzp*1.5*dx
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
      p_i = solnData(TEMP_VAR,ielem(i,IAXIS,presflag+1), &
                              ielem(i,JAXIS,presflag+1), &
                              ielem(i,KAXIS,presflag+1));
      zpres2 = zpres2 + phile(i,1)*p_i
  enddo
  ! Release Pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  ! Substract hydrostatic pressure to total pressure in point 2:
  zpres2 = zpres2

  !! Point 3:
  xbe3(IAXIS) = xp
  xbe3(JAXIS) = yp
#if NDIM == 3
  xbe3(KAXIS) = zp
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
  zpres3 = zpres3 


  particleData(TL_PART_PROP)   = zpres3
  particleData(HFLX_PART_PROP) = (zpres3 - zpres2)/(1.5*dx*zpres3)

  return

end subroutine ib_distributedFluxes


