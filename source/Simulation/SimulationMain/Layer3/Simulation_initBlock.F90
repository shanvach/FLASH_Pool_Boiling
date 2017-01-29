!!****if* source/Simulation/SimulationMain/Layer3/Simulation_initBlock
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
!! This problem is that of an intense
!! laser hitting a 3-layer target of decreasing density. The
!! three materials in the target are copper, CH plastic, and 
!! C-foam. The first interface, between the copper and plastic
!! has a sinusoidal perturbation, and this perturbation is
!! imprinted on the shock that passes through the interface. As
!! the shock propagates toward the second interface, the perturbation
!! oscillates. When the shock hits the second interface, it
!! is almost planar, but the perturbation is still imprinted on
!! to the interface. The interfaces are subject to fluid instabitlites
!! resulting from the passage of the shock, and instability growth
!! (fingers of copper and bubbles of C-foam) may is observed.
!! 
!! Reference: Kane et al., Phys. Rev. E, 63 055401(R) 2001
!!
!!
!!
!! ARGUMENTS
!!      blockID      -       the current block number to be filled
!!      myPE         -       the current processor number 
!!
!!
!!***
subroutine Simulation_initBlock(blockId)
  use Simulation_data, ONLY : sim_pi,&
                              sim_xmin, sim_xmax, sim_ymin, sim_ymax, &
                              sim_small,sim_xzn1d,sim_1dModel
  use Eos_interface, ONLY : Eos

  use Simulation_data, ONLY :sim_dens = 1,sim_velx = 2,sim_vely = 3,&
                             sim_velz = 4,sim_pres = 5,sim_ener = 6,&
                             sim_temp = 7,sim_gamc = 8,sim_game = 9,&
                             sim_enuc = 10,sim_gpot = 11,sim_speciesBegin = 12
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  integer, intent(IN) :: blockID, myPE

  integer :: ierr
  integer, parameter ::  ifail = -1
  
  real :: dxinit, dx_fine
  real :: xsh, xin1,xin2,xin3
  integer :: izone,izonecu,izonech

  real :: abar, zbar, ytemp
  real :: ylength
  
  integer :: status
  
  
  
  !     the perturbation parameters
  real :: amplitude,lambda
  real :: xpert
  real :: pcu,eicu,gamcu,gamecu,tcu
  real :: pch,eich,gamch,gamech,tch
  real :: rjunk1,rjunk2
  real :: xshift,xfudge
  integer, dimension(LOW:HIGH,MDIM):: blkLimits,blkLimitsGC
  integer :: sizeX, sizeY, sizeZ
  real, allocatable, DIMENSION(:) :: rho, p, ei, ek, tmp, gamc,   & 
          u, v, w     !storage for 0's to be compatible with dBasePutData 
  real, allocatable, DIMENSION(:,:) :: xn 

  real, allocatable, DIMENSION(:) :: xLeft,yLeft,zLeft,&
                                     xCenter,yCenter,zCenter, &
                                     xRight, yRight, zRight
  real :: tmp_velx, tmp_vely, tmp_velz ! temporary storage for velocities
  real, dimenstion(:,:,:,:)::sVec    
!!$      real, DIMENSION(nvar, iLo_gc:iHi_gc, & 
!!$     &     jLo_gc:jHi_gc,kLo_gc:kHi_gc) ::  sVec
  integer :: istat
  real, dimension(EOS_NUM) :: eosData
  logical, dimension(EOS_VARS+1:EOS_NUM) :: mask
  real, dimension(NSPECIES) :: eosMassFr
  integer :: vecLen=1

!! allocate space for the vectors
  call Grid_getBlkLimits(blockID, blkLimits, blkLimitsGC)
  sizeX=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
  allocate(rho(sizeX),stat=istat)
  allocate(p(sizeX),stat=istat)
  allocate(ei(sizeX),stat=istat)
  allocate(ek(sizeX),stat=istat)
  allocate(tmp(sizeX),stat=istat)
  allocate(gamc(sizeX),stat=istat)
  allocate(u(sizeX),stat=istat)
  allocate(v(sizeX),stat=istat)
  allocate(w(sizeX),stat=istat)
  allocate(xLeft(sizeX),stat=istat)
  allocate(xCenter(sizeX),stat=istat)
  allocate(xRight(sizeX),stat=istat)
  allocate(yLeft(sizeY),stat=istat)
  allocate(yCenter(sizeY),stat=istat)
  allocate(yRight(sizeY),stat=istat)
  allocate(zLeft(sizeZ),stat=istat)
  allocate(zCenter(sizeZ),stat=istat)
  allocate(zRight(sizeZ),stat=istat)
  allocate(xn(sizeX,NSPECIES),stat=istat)
  allocate(sVec(NUNK_VARS,sizeX,sizeY,sizeZ),stat=istat)
! set the shock and unperturbed interface x positions

  xshift = 1.5e-2
!     shift to fix Jave's initial conditions
!     2d      xfudge = 1.5e-3
!     1d      xfudge = 0.0e0
  xfudge = 1.5e-3
!     things measured as read from the 1-d input
  xsh  = 0.0221e0           !shock
  xin2 = 0.0227e0           !Cu-CH
  xin3 = 0.0392e0           !CH-Cfoam
!     dxinit = 0.01E-01
  dxinit = (xin2-xsh)/2.0e0
      
! set the perturbation parameters
  lambda = 0.02e0
!     2d      amplitude = 0.0015e0
!     1d      amplitude = 0.0000e0
  amplitude = 0.0015e0

     
  dx_fine = (sim_xmax - sim_xmin)/real(N1D)

  ylength = sim_ymax - sim_ymin
     
! pick a zone that is slightly below the Cu-CH and reference it as the cu zone
  izonecu = int(((xin2 - dxinit) - sim_xmin)/dx_fine) + 1

! set the mass fraction (should be pure cu)
  do n = 1, NSPECIES
     xn(1,n) = sim_1dModel(izonecu,sim_speciesBegin-1+n)
  enddo
  mask = .false.
  eos_data(EOS_DENS)=sim_1dModel(izonecu,sim_dens)
  eos_data(EOS_EINT)=sim_1dModel(izonecu,sim_ener)
  call Eos(MODE_DENS_EI,vecLen,eosData,eosMassFr,mask)
  tcu=eosData(EOS_TEMP)
  gammac=eosData(EOS_GAMC)
  pcu=eosData(EOS_PRES)

! call the eos to set tcu, the temperature in the copper region
  call Eos (model_1d(izonecu,sim_dens), tcu, pcu,  & 
     &          model_1d(izonecu,sim_ener), xn(i,:),entropy,  abar, zbar, & 
     &          dpt, dpd, det, ded,dst,dsd,   c_v, c_p, gammac, pel, ne, eta, & 
     &          2)

! pick a zone that is slightly above the 2nd interface and 
! reference it as the ch zone
      izonech = int(((xin2 + 3.0e0*dxinit) - sim_xmin)/dx_fine) + 1

! set the gamma for ch
      call get_fluid_property("ch","adiabatic index",gamch,status)
     
! set the mass fraction (should be pure ch)
      do n = 1, NSPECIES
         xn(i,n) = model_1d(izonech,sim_speciesBegin-1+n)
      enddo

! call the eos to set tch, the temperature in the ch region
      call Eos (model_1d(izonech,sim_dens), tch, pch,  & 
     &          model_1d(izonech,sim_ener), xn(i,:),entropy,  abar, zbar, & 
     &          dpt, dpd, det, ded,dst,dsd,   c_v, c_p, gammac, pel, ne, eta, & 
     &          2)

! initially no y or z velocity
      v(:) = 0.0e0
      w(:) = 0.0e0
!
! set coordinages
!
      y = 0.
      yl = 0.
      yr = 0.
      z = 0.
      zl = 0.
      zr = 0.
      if (ndim > 2) then
         call dBaseGetCoords(iznr, iZcoord, block_no, zr)
         call dBaseGetCoords(izn , iZcoord, block_no, z )
         call dBaseGetCoords(iznl, iZcoord, block_no, zl)
      endif
      if (ndim > 1) then
         call dBaseGetCoords(iznr, iYcoord, block_no, yr)
         call dBaseGetCoords(izn , iYcoord, block_no, y )
         call dBaseGetCoords(iznl, iYcoord, block_no, yl)
      endif
      call dBaseGetCoords(iznr, iXcoord, block_no, xr)
      call dBaseGetCoords(izn , iXcoord, block_no, x )
      call dBaseGetCoords(iznl, iXcoord, block_no, xl)

! now fill the master arrays
      do k = nguard*k3d+1, nguard*k3d+nzb
         do j = nguard*k2d+1, nguard*k2d+nyb
            do i = nguard+1, nguard+nxb

! find the x zone in the 1-d array
               izone = int((x(i) - sim_xmin)/dx_fine) + 1

! fill the scratch arrays and 

               tmp(i) = model_1d(izone,sim_temp)
               rho(i) = model_1d(izone,sim_dens)
               ei(i)  = model_1d(izone,sim_ener)
               p(i)   = model_1d(izone,sim_pres)
               u(i)   = model_1d(izone,sim_velx)

               do n = 1, ionmax
                  xn(i,n) = model_1d(izone,sim_speciesBegin-1+n)
               enddo
! call the eos to get it consistent    

      call Eos (rho(i), tmp(i), p(i), & 
     &          ei(i), xn(i,:),entropy,  abar, zbar, & 
     &          dpt, dpd, det, ded,dst,dsd,   c_v, c_p, gammac, pel, ne, eta, & 
     &          2)

! now store things according to where they are relative to the
! to the perturbation.

               if ((x(i).le.(xsh+dxinit)).or.(x(i).ge.xin3)) then
! below shock + eps or above interface + amplitude of the perturbation
! just stuff the arrays
                  call dBasePutData(idens, "Point", i, j, k,  & 
     &                 block_no, rho(i), 1)
                  call dBasePutData(ivelx, "Point", i, j, k,  & 
     &                 block_no, u(i), 1)
                  call dBasePutData(ipres, "Point", i, j, k,  & 
     &                 block_no, p(i), 1)
                  call dBasePutData(iener, "Point", i, j, k, block_no,   & 
     &                 model_1d(izone,sim_ener) +  & 
     &                 .5*u(i)**2, 1 )
                  call dBasePutData(itemp, "Point", i, j, k, block_no,   & 
     &                 tmp(i), 1)
                  call dBasePutData(igamc, "Point", i, j, k, block_no,   & 
     &                 gammac, 1)
                  call dBasePutData(igame, "Point", i, j, k, block_no,   & 
     &                 p(i)/(ei(i)*rho(i)) + 1.0, 1 )

                  if (x(i).gt.xin3) then
                     call dBasePutData(idbcu, "Point", & 
     &                     i, j, k, block_no, small, 1)
                     call dBasePutData(idbch, "Point", & 
     &                     i, j, k, block_no, small, 1)
                     call dBasePutData(idbcf, "Point", & 
     &                     i, j, k, block_no, 1.0e0 - small, 1)
                  else
                     call dBasePutData(idbcu, "Point", & 
     &                     i, j, k, block_no, 1.0e0 - small, 1)
                     call dBasePutData(idbch, "Point", & 
     &                     i, j, k, block_no, small, 1)
                     call dBasePutData(idbcf, "Point", & 
     &                     i, j, k, block_no, small, 1)
                  endif
               else
!  set the location of the interface
                  xpert = xin2 + xfudge + amplitude*sin(2.0e0*sim_pi & 
                       &                     *(y(j)-sim_ymin)/lambda - (3.0e0*sim_pi)/2.0e0)
!  if less than the perturbation, then cu
                  if(x(i).lt.xpert) then
                     call dBasePutData(idens, "Point", i, j, k, & 
     &                    block_no, model_1d(izonecu,sim_dens), 1)
                     call dBasePutData(ivelx, "Point", i, j, k, & 
     &                    block_no, model_1d(izonecu,sim_velx), 1)
                     call dBasePutData(ipres, "Point", i, j, k, & 
     &                    block_no, pcu, 1)
                     call dBasePutData(iener, "Point", i, j, k,  & 
     &                    block_no, model_1d(izonecu,sim_ener) + & 
     &                    .5*model_1d(izonecu,sim_velx)**2, 1 )
                     call dBasePutData(itemp, "Point", i, j, k, & 
     &                    block_no, tcu, 1)
                     call dBasePutData(igamc, "Point", i, j, k, & 
     &                    block_no, sim_gamcu, 1)
                     call dBasePutData(igame, "Point", i, j, k, & 
     &                    block_no, & 
     &                    model_1d(izonecu,sim_pres)/ & 
     &                    (model_1d(izonecu,sim_ener)* & 
     &                    model_1d(izonecu,sim_dens)) + 1.0e0, 1 )

                     if (x(i).lt.xin3) then
                        call dBasePutData(idbcu, "Point", & 
     &                         i, j, k, block_no, 1.0e0 - small, 1)
                        call dBasePutData(idbch, "Point", & 
     &                         i, j, k, block_no, small, 1)
                        call dBasePutData(idbcf, "Point", & 
     &                         i, j, k, block_no, small, 1)
                     endif

!  else greater than the perturbation, then ch
                  else
                     call dBasePutData(idens, "Point", i, j, k, & 
     &                    block_no, model_1d(izonech,sim_dens), 1)
                     call dBasePutData(ivelx, "Point", i, j, k, & 
     &                    block_no, model_1d(izonech,sim_velx), 1)
                     call dBasePutData(ipres, "Point", i, j, k, & 
     &                    block_no, pch, 1)
                     call dBasePutData(iener, "Point", i, j, k, & 
     &                    block_no, model_1d(izonech,sim_ener) + & 
     &                    .5*model_1d(izonech,sim_velx)**2, 1 )
                     call dBasePutData(itemp, "Point", i, j, k, & 
     &                    block_no, tch, 1)
                     call dBasePutData(igamc, "Point", i, j, k, & 
     &                    block_no,gamch, 1)
                     call dBasePutData(igame, "Point", i, j, k,  & 
     &                    block_no, & 
     &                    model_1d(izonech,sim_pres)/ & 
     &                    (model_1d(izonech,sim_ener)* & 
     &                    model_1d(izonech,sim_dens)) + 1.0e0, 1 )

                     if (x(i).lt.xin3) then
                       call dBasePutData(idbcu, "Point", & 
     &                        i, j, k, block_no, small, 1)
                       call dBasePutData(idbch, "Point", & 
     &                        i, j, k, block_no, 1.0e0 - small, 1)
                       call dBasePutData(idbcf, "Point", & 
     &                        i, j, k, block_no, small, 1)
                     endif
                  endif
               endif
!
            enddo

!     Dbase can be used to store a single point or 1/2/3d arrays of data. 
!     I chose to do y,z-velocities as 1-d arrays, though I could have 
!     chosen to store them point by point as well. 
            call dBasePutData("vely", "xVector", j, k, block_no, v, 1 ) 
            call dBasePutData("velz", "xVector", j, k, block_no, w, 1 )
            
         enddo
      enddo
!     
      call dBaseGetData("AllVariables", block_no, sVec, 1)
!     
!     
!     as a test, write out a strip
      if(myPE.eq.0) then
        open(unit=89,file='t.out',status='unknown',position & 
     &        ='append')
        do i =  nguard+1, nguard+nxb
            write(89,102) x(i) & 
     &     ,sVec(idens,i,nguard*k2d+1,nguard*k3d+1) & 
     &     ,sVec(ivelx,i,nguard*k2d+1,nguard*k3d+1) & 
     &     ,sVec(ively,i,nguard*k2d+1,nguard*k3d+1) & 
     &     ,sVec(ivelz,i,nguard*k2d+1,nguard*k3d+1) & 
     &     ,sVec(ipres,i,nguard*k2d+1,nguard*k3d+1) & 
     &     ,sVec(iener,i,nguard*k2d+1,nguard*k3d+1) & 
     &     - 0.5e0*(sVec(ivelx,i,nguard*k2d+1,nguard*k3d+1)**2) & 
     &     ,sVec(igamc,i,nguard*k2d+1,nguard*k3d+1) & 
     &     ,sVec(igame,i,nguard*k2d+1,nguard*k3d+1) & 
     &     ,sVec(idbcu,i,nguard*k2d+1,nguard*k3d+1) & 
     &     ,sVec(idbch,i,nguard*k2d+1,nguard*k3d+1) & 
     &     ,sVec(idbcf,i,nguard*k2d+1,nguard*k3d+1)
        enddo
        write(89,*)
        close(89)
      endif

 102  format(1E14.8,2x,1E14.8,2x,1E14.8,2x,1E14.8,2x,1E14.8,2x,1E14.8,2x, & 
     &     1E14.8,2x,1E14.8,2x,1E14.8,2x,1E14.8,2x,1E14.8,2x,1E14.8,2x, & 
     &     1E14.8,2x,1E14.8,2x,1E14.8, & 
     &     2x,1E14.8,2x,1E14.8,2x,1E14.8)


  deallocate(rho)
  deallocate(p)
  deallocate(ei)
  deallocate(ek)
  deallocate(tmp)
  deallocate(gamc)
  deallocate(u)
  deallocate(v)
  deallocate(w)
  deallocate(xLeft)
  deallocate(xCenter)
  deallocate(xRight)
  deallocate(yLeft)
  deallocate(yCenter)
  deallocate(yRight)
  deallocate(zLeft)
  deallocate(zCenter)
  deallocate(zRight)
  deallocate(xn)
  deallocate(sVec)
  return
end subroutine Simulation_initBlock
