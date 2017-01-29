!!****if* source/Simulation/SimulationMain/ShockCyl/Simulation_initBlock
!!
!! NAME
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!  Simulation_initBlock(  integer(in) :: blockID,
!!                         
!!
!! DESCRIPTION   
!!     Initialize fluid properties (density, pressure, velocity, etc.) 
!!          in a single block 
!!          for the 1 cylinder Shock Cylinder simulation in 2D or 3D
!!
!! ARGUMENTS
!!      blockID    integer:in      the current block number to be filled
!!      myPE        integer:in      the current processor number 
!!
!!
!!***


subroutine Simulation_initBlock(blockId)

!============================================================================

 
  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_getDeltas, Grid_putPointData, &
    Grid_putRowData
  use Eos_interface, ONLY : Eos
!!  use physicaldata, ONLY : unk

  implicit none 



#include "constants.h"
#include "Flash.h"
#include "Eos.h"


  integer, intent(IN) :: blockID, myPE

  ! declare the eos arguments
  real :: gammac,  pres
  real :: eint, entropy, gamsf6, t, dst, dsd, gamma, xxne, xalfa 
  real :: gamair, abar, zbar, dist



  integer               :: i, j, k, n, status,numXn,istat
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord,dx,dy,dz
  real, allocatable,dimension(:)    :: rho, p, ei, ek, vx, vt, ge, gc, e
  real, allocatable,dimension(:)    :: vy, vz, game, gamc

  real, dimension(NSPECIES)   :: xnshock, xnamb

  real            :: gammam1, gam1inv


  real       :: frac,eps
  integer    :: ix,iy
  real       :: xx, yy, zz
  real       :: gp1, gm1
  real       :: vul, vur, vll, vlr, termx, termy
  real       :: vxu, vxl, vcnew, denom 
  real       :: g1, g2, delp, vmin, func, dr
  real       :: mw_mix

  real       :: radius
  real       :: x_sf6, rvel, zvel
  real       :: rz_pert_fact
  real       ::  x_loc,  y_loc,  z_loc, r_loc
  real       :: dx_loc, dy_loc, dz_loc
  real       :: x_sf6_sum, rvel_sum, zvel_sum, deni
  integer    :: ii, jj, kk, vecLen =1

  real, dimension(NSPECIES) :: massFrac

  integer    :: iylow, iyhigh, ixlow, ixhigh, l, m
  integer :: szx, szy, szz
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: pos
  real, dimension(MDIM) :: delta
  integer :: dataSize
  !!DEV: kda verify this.  temporary for now




  
  !----------------------------------------------------------------------------
  
  ! get integer index sizes of the block
  ! blkLimits holds the interior block limits, blkLimitsGC holds entire block 
  ! indices including guardcells
  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
  szx=blkLimitsGC(HIGH,IAXIS)
  szy=blkLimitsGC(HIGH,JAXIS)
  szz=blkLimitsGC(HIGH,KAXIS)

  numXn = max(szx,szy)
  numXn = max(numXn,szz)

  allocate(rho(numXn),stat=istat)
  allocate(p(numXn),stat=istat)
  allocate(ei(numXn),stat=istat)
  allocate(ek(numXn),stat=istat)
  allocate(vx(numXn),stat=istat)
  allocate(vy(numXn),stat=istat)
  allocate(vz(numXn),stat=istat)
  allocate(vt(numXn),stat=istat)
  allocate(ge(numXn),stat=istat)
  allocate(gc(numXn),stat=istat)
  allocate(e(numXn),stat=istat)
  allocate(game(numXn),stat=istat)
  allocate(gamc(numXn),stat=istat)


  allocate(xCoord(szx),stat=istat)
  allocate(dx(szx),stat=istat)
  allocate(yCoord(szy),stat=istat)
  allocate(dy(szy),stat=istat)
  allocate(zCoord(szz),stat=istat)
  allocate(dz(szz),stat=istat)


! test results without guard cells
!  gcell = .false.

  call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell,xCoord,szx)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell,yCoord,szy)
  call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell,zCoord,szz)
  gcell = .true.

  call Grid_getDeltas(blockID,delta)
  dx = delta(IAXIS)
  dy = delta(JAXIS)
  dz = delta(KAXIS)
!  print *,'x=',xCoord
!  print *,'dx=',dx
!  print *,'y=',yCoord
!  print *,'dy=',dy


!----------------------------------------------------------------------------

! Initialize the flowfield.
! First, each cell is initialized as either pre-shock or post-shock air;
!   then the cylinder conditions are superimposed, if appropriate.

! determine the post-shock conditions fro the ambient (amb) conditions. 
! 'ps' refers to post-shock
! 'p' stands for pressure, 'rho' is density and 'v' is velocity


  mode = MODE_DENS_PRES   ! requires density and pressure as input, returns others
  eos_arr(EOS_PRES)=p_amb
  eos_arr(EOS_DENS)=rho_amb

  massFrac(species_air) = 1.e0 - sim_smallSF6
  massFrac(species_sf6) =        sim_smallSF6

  call Eos(mode, vecLen, eos_arr, massFrac)

  t = eos_arr(EOS_TEMP)
  eint = eos_arr(EOS_EINT)
  gp1    = eos_arr(EOS_GAMC) + 1.e0
  gm1    = eos_arr(EOS_GAMC) - 1.e0
  g1     = gp1/gm1
  g2     = 2.e0*eos_arr(EOS_GAMC)/gp1
  delp   = 1.e0 + g2*(mach**2 - 1.e0)


  p_ps   = p_amb * delp
  rho_ps = rho_amb * (1.e0 + g1*delp)/(g1+delp)

  vx_ps  = (sqrt(eos_arr(EOS_GAMC)*p_amb/rho_amb)/eos_arr(EOS_GAMC))*(delp -1.e0)
  vx_ps  = vx_ps * sqrt(g2/(delp + 1.e0/g1))

  ! move some initializations outside of loop -- should be moved to Simulation_init eventually
  if (use_rz_sim_data) then 
      rz_pert_fact = 2.0*pi/rz_pert_zlen 
  end if 
  


  ! Loop over all cells (including guardcells) in the block and initialize the variables.
  
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     zz = zCoord(k)
     
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        yy = yCoord(j)
        
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           xx = xCoord(i)
           
           ! to be extra safe, initialize mass fractions for every
           ! cell.  some problems occur if we don't because we're 
           ! passing over guardcells, too
           massFrac(species_air) = 1.e0 - sim_smallSF6
           massFrac(species_sf6) =        sim_smallSF6
           
           ! Initialize cells to the left of the shock.
           
           if ( xCoord(i) < sim_xShock ) then ! post-shock region
              
              rho(i)   = rho_ps
              p(i)     = p_ps
              vx(i)    = vx_amb + vx_ps
              vy(i)    = 0.e0
              vz(i)    = 0.e0
              
           else ! ambient medium
              
              rho(i)   = rho_amb
              p(i)     = p_amb
              vx (i)   = vx_amb
              vy(i)    = 0.e0
              vz(i)    = 0.e0

           end if

           ! SF6 cylinder

           if ( sim_useRawData ) then

              
              ixlow  = sim_rawNumPixelsX * ((xx - sizex(1))/(sizex(sim_rawNumPixelsX)-sizex(1)))
              ixhigh = ixlow + 1

              iylow  = sim_rawNumPixelsY * ((yy - sizey(1))/(sizey(sim_rawNumPixelsY)-sizey(1)))
              iyhigh = iylow + 1

              if ( ixlow.ge.1 .and. ixhigh.le.sim_rawNumPixelsX .and. &
                   iylow.ge.1 .and. iyhigh.le.sim_rawNumPixelsY ) then

                 vul = vconc(ixlow,iyhigh)
                 vur = vconc(ixhigh,iyhigh)
                 vll = vconc(ixlow,iylow)
                 vlr = vconc(ixhigh,iylow)

                 termx = (xx - sizex(ixlow))/(sizex(ixhigh) - sizex(ixlow))
                 termy = (yy - sizey(iylow))/(sizey(iyhigh) - sizey(iylow))

                 vxl = vll + (vlr - vll) * termx
                 vxu = vul + (vur - vul) * termx

                 vcnew = vxl + (vxu - vxl) * termy

                 vmin = min(vul,vur,vll,vlr)
                 vmax = max(vul,vur,vll,vlr)

                 vcnew = max(vmin,sim_smallSF6,min(vmax, vcnew ))

!                vcnew is the mole fraction of sf6.
!                now get mass fractions:
                 mw_mix = vcnew*mw_sf6 + (1.e0-vcnew)*mw_air

                 massFrac(species_sf6) = vcnew*mw_sf6/mw_mix
                 massFrac(species_air) = (1.e0 - vcnew)*mw_air/mw_mix

                 rho(i) = c_initial*mw_mix
                 p(i)   = p_amb
                 vx(i)  = vx_amb
                 vy(i)  = 0.e0
                 vz(i)  = massFrac(species_sf6)*vz_fact
                 
              end if
              
           elseif (use_rz_sim_data) then  ! data from rz simulation data
              

              dz_loc = dz(k)/float(rz_subintNZ)
              dy_loc = dy(j)/float(rz_subintNY)
              dx_loc = dx(i)/float(rz_subintNX)
              
              radius   = sqrt((xx-xctr)**2 + (yy-yctr)**2)
!              print *,'i,xx,radius=',i,xx,radius
!              print *,'j,yy       =',j,yy
!              print *,'xx,yy,radius = ',xx,yy,radius
              if ( (  radius       <= rz_rmax-max(dx(i), dy(j))) .and. &
                   (  xx           > sim_xShock )                       )  then

                 !                Perturbations applied to evaluation location of radius
                 !                Only applied for z-direction symmetry in 3d, in 
                 !                  which case the z-location is rz_zplane
                 
                 
                 !                Numerical integration (average) over the cell to
                 !                  get x_sf6, rvel, and if 3d, zvel.
                 x_sf6_sum = 0.0
                 rvel_sum  = 0.0
                 zvel_sum  = 0.0
                 deni = 1.0/float(rz_subintNX * rz_subintNY * rz_subintNZ)
!           print *,'z_loc= ',z_loc,' zz= ',zz,' dz =',dz(k),' dz_loc= ',dz_loc, ' k,kk= ',k,kk      
                 do ii = 1, rz_subintNX
                    do jj = 1, rz_subintNY
                       
                       if (NDIM == 3) then
                          do kk = 1, rz_subintNZ
                             x_loc = xx - 0.5*dx(i) + (float(ii-1) + 0.5)*dx_loc
                             y_loc = yy - 0.5*dy(j) + (float(jj-1) + 0.5)*dy_loc
                             z_loc = zz - 0.5*dz(k) + (float(kk-1) + 0.5)*dz_loc
                             r_loc = sqrt( (x_loc - xctr)**2 + (y_loc - yctr)**2)
           
                             
                             if(rz_3d_use_sym) then
                                r_loc = abs(r_loc + rz_pert_amp*             &
                                     sin(z_loc*rz_pert_fact))
                                r_loc = min(r_loc,rz_rmax)
                                z_loc = rz_zplane 
                             endif
                             
                             call sim_rzInitialConditions(r_loc, z_loc, x_sf6, rvel, zvel)
                             
                             x_sf6_sum = x_sf6_sum + x_sf6
                             rvel_sum  = rvel_sum + rvel
                             zvel_sum  = zvel_sum + zvel
                          enddo  ! kk
                       else  ! NDIM == 2
                          x_loc = xx - 0.5*dx(i) + (float(ii-1) + 0.5)*dx_loc
                          y_loc = yy - 0.5*dy(j) + (float(jj-1) + 0.5)*dy_loc
                          r_loc = sqrt( (x_loc - xctr)**2 + (y_loc - yctr)**2)
                          z_loc = rz_zplane
                          
                          call sim_rzInitialConditions(r_loc, z_loc, x_sf6, rvel, zvel)
                          
                          x_sf6_sum = x_sf6_sum + x_sf6
                          rvel_sum  = rvel_sum + rvel
                          !                      zvel_sum remains 0.0, as initialized
                       endif ! NDIM
                       
                    enddo  ! jj
                 enddo   ! ii
                 x_sf6= x_sf6_sum*deni
                 rvel = rvel_sum*deni
                 zvel = zvel_sum*deni
                 
                 !                Limit sf6 mass fraction
                 massFrac(species_sf6) = max(sim_smallSF6, min(1.0e0-sim_smallSF6, x_sf6))
                 massFrac(species_air) = 1.0e0 - massFrac(species_sf6)
                 !                Set pressure to ambient pressure - in rz_sim,
                 !                  the pressure is dominated by hydrostatic.
                 p(i) = p_amb
                 
                 vx(i)  = vx_amb + rvel*(xx-xctr)/radius
                 vy(i)  =          rvel*(yy-yctr)/radius
                 vz(i)  = zvel
                 
                 mw_mix = 1.0/( massFrac(species_sf6)/mw_sf6 + massFrac(species_air)/mw_air )
                 rho(i) = c_initial*mw_mix
              end if

              

              
           elseif (sim_useRadialFit)  then ! data from radial fit formula
              
          
              radius   = sqrt((xx-xctr)**2 + (yy-yctr)**2)
              if ( (radius <= sim_radialFitRadius)  .and.     &
                   (xx     > sim_xShock   )         )  then
                 call sim_radialFit(radius,maxconc,vcnew)
                 !                vcnew is the mole fraction of sf6.
                 !                now get mass fractions:
                 mw_mix = vcnew*mw_sf6 + (1.e0-vcnew)*mw_air
                 massFrac(species_sf6) = vcnew*mw_sf6/mw_mix
                 massFrac(species_air) = (1.e0 - vcnew)*mw_air/mw_mix
                 
                 rho(i) = c_initial*mw_mix
                 p(i)   = p_amb
                 vx(i)  = vx_amb
                 vy(i)  = 0.e0
                 vz(i)  = massFrac(species_sf6)*vz_fact
              end if
              
           end if  ! raw data, rz_sim_data, or radial fit
           
           
           !-------------------------------------------------------------------
           ! Calculate temp and internal energy from equation of state
           !-------------------------------------------------------------------

           eos_arr(EOS_DENS) = rho(i)
           eos_arr(EOS_PRES) = p(i)
           call Eos(mode, vecLen, eos_arr, massFrac)

           !print *, "massFrac = ", massFrac

           eint=eos_arr(EOS_EINT)
           gamma = eos_arr(EOS_GAMC)
           t = eos_arr(EOS_TEMP)

           e(i) = eint + 0.5e0*(vx(i)**2 + vy(i)**2 + vz(i)**2)

           game(i) = p(i)/(eint*rho(i)) + 1.e0
           gamc(i) = gamma

           !-------------------------------------------------------------------
           ! finally, fill the solution array
           !-------------------------------------------------------------------
           pos(IAXIS) = i
           pos(JAXIS) = j
           pos(KAXIS) = k

           do n=SPECIES_BEGIN,SPECIES_END

              call Grid_putPointData(blockID, CENTER, n, EXTERIOR, pos, massFrac(n-NPROP_VARS))

           enddo
           
        enddo  ! sweep over x

        pos(IAXIS) = 1
        pos(JAXIS) = j
        pos(KAXIS) = k
        dataSize=blkLimitsGC(HIGH,IAXIS)

  
        call Grid_putRowData(blockID, CENTER, DENS_VAR, EXTERIOR, IAXIS, pos, rho, dataSize)
        call Grid_putRowData(blockID, CENTER, PRES_VAR, EXTERIOR, IAXIS, pos, p, dataSize)
        call Grid_putRowData(blockID, CENTER, ENER_VAR, EXTERIOR, IAXIS, pos, e, dataSize)    
        call Grid_putRowData(blockID, CENTER, GAME_VAR, EXTERIOR, IAXIS, pos, gamc, dataSize)
        call Grid_putRowData(blockID, CENTER, GAMC_VAR, EXTERIOR, IAXIS, pos, game, dataSize)
        call Grid_putRowData(blockID, CENTER, VELX_VAR, EXTERIOR, IAXIS, pos, vx, dataSize)
        call Grid_putRowData(blockID, CENTER, VELY_VAR, EXTERIOR, IAXIS, pos, vy, dataSize)
        call Grid_putRowData(blockID, CENTER, VELZ_VAR, EXTERIOR, IAXIS, pos, vz, dataSize)





     enddo     ! sweep over y
  enddo        ! sweep over z



  !============================================================================
  deallocate(rho)
  deallocate(p)
  deallocate(ei)
  deallocate(ek)
  deallocate(vx)
  deallocate(vy)
  deallocate(vz)
  deallocate(vt)
  deallocate(ge)
  deallocate(gc)
  deallocate(e)
  deallocate(game)
  deallocate(gamc)


  deallocate(xCoord)
  deallocate(dx)
  deallocate(yCoord)
  deallocate(dy)
  deallocate(zCoord)
  deallocate(dz)

  return
end subroutine Simulation_initBlock
