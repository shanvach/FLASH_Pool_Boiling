!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkAccelSinksOnGas
!!
!! NAME
!!
!!  Particles_sinkAccelSinksOnGas
!!
!! SYNOPSIS
!!
!!  call Particles_sinkAccelSinksOnGas(integer,intent(IN)  :: blockcount,
!!                                     integer,dimension(blockCount),intent(IN)  :: blocklist)
!!
!! DESCRIPTION
!!
!!  Computes SGAX,SGAY,SGAY unk vars from sink particles.
!!
!! ARGUMENTS
!!
!!   blockcount : 
!!
!!   blocklist : 
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!
!!***

! used to be SinksOnGasAccel
subroutine Particles_sinkAccelSinksOnGas(blockCount,blockList)

!==============================================================================

 use Particles_sinkData
 use Driver_data, ONLY : dr_globalMe
 use Driver_interface, ONLY : Driver_abortFlash
 use Grid_interface, ONLY : Grid_getCellCoords, Grid_getBlkIndexLimits,  & 
     Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkPhysicalSize
 use RuntimeParameters_interface, ONLY : RuntimeParameters_get
 use Cosmology_interface, ONLY : Cosmology_getRedshift
 use PhysicalConstants_interface, ONLY : PhysicalConstants_get

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

 implicit none

 integer,intent(IN) :: blockCount
 integer,dimension(blockCount),intent(IN) :: blockList

 real, POINTER, DIMENSION(:,:,:,:) :: solnData

 character(len=80), save :: softening_type_gas, grav_boundary_type

 logical, save :: first_call = .true.
 real, save    :: newton, maxradius_pbc
 integer, save :: softeningtype
 real, save    :: softening_radius
 real          :: slope, hinv, h2inv
 real, save    :: xmin, xmax, ymin, ymax, zmin, zmax, Lx, Ly, Lz

 real          :: radius, x2, y2, z2, prefactor, q, kernelvalue, r3
 real          :: paccx, paccy, paccz, mass
 integer       :: nxbBlock, nybBlock, nzbBlock
 real          :: blockSize(MDIM)
 real          :: dx, dy, dz, dx_inside, dy_inside, dz_inside
 integer       :: ii,jj,kk, p, nx, ny, nz

 real, dimension(:), allocatable :: x, y, z
 integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
 integer       :: size_x, size_y, size_z
 integer       :: lb
 real          :: redshift, softening_radius_comoving, oneplusz3

 !==============================================================================

 if (first_call) then

   call PhysicalConstants_get("Newton", newton)

   if (UseSinkParticles) then
      call RuntimeParameters_get("sink_softening_radius", softening_radius)
      call RuntimeParameters_get("sink_softening_type_gas", softening_type_gas)
      select case (softening_type_gas)
         case ("spline")
            softeningtype = 1
         case ("linear")
            softeningtype = 2
         case default
            softening_type_gas = "linear"
            softeningtype = 2
            if (dr_globalMe .eq. MASTER_PE) print*, 'invalid sink_softening_type_gas specified. using default-> ', &
                                              & trim(softening_type_gas)
      end select
      if (dr_globalMe .eq. MASTER_PE) print*, 'Particles_sinkAccelSinksOnGas:: sink_softening_type_gas = ', &
                                        & trim(softening_type_gas)

      call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)

      if ((grav_boundary_type .ne. "isolated") .and. (grav_boundary_type .ne. "periodic")) then
        call Driver_abortFlash('Sink particles only work with isolated or periodic gravity boundaries.')
      endif

      call RuntimeParameters_get('xmin', xmin)
      call RuntimeParameters_get('xmax', xmax)
      call RuntimeParameters_get('ymin', ymin)
      call RuntimeParameters_get('ymax', ymax)
      call RuntimeParameters_get('zmin', zmin)
      call RuntimeParameters_get('zmax', zmax)

      Lx = xmax-xmin
      Ly = ymax-ymin
      Lz = zmax-zmin
      if (grav_boundary_type .eq. "periodic") then
         maxradius_pbc = real(nrep_pbc)*min(Lx,Ly,Lz)
         if (dr_globalMe .eq. MASTER_PE) print*, 'Particles_sinkAccelSinksOnGas:: maxradius_pbc = ', maxradius_pbc
      endif
   else
      softening_radius = 0.e0
      slope  = 0.e0
   endif

   first_call = .false.
 endif

 !==============================================================================

 if (localnpf .EQ. 0) return

 call Cosmology_getRedshift(redshift)
 softening_radius_comoving = softening_radius * (1.0 + redshift)
 hinv  = 2.0/softening_radius_comoving !!! makes sure that only for r < r_soft actual softening occurs
 h2inv = hinv**2
 slope = 1.0/softening_radius_comoving**3
 oneplusz3 = (1.0 + redshift)**3.0
 prefactor = -newton*oneplusz3

 if (grav_boundary_type .eq. "periodic") then

   ! Loop through blocks
   do lb = 1, blockCount

      call Grid_getBlkPtr(blockList(lb),solnData)

      call Grid_getBlkIndexLimits(lb, blkLimits, blkLimitsGC)
      call Grid_getBlkPhysicalSize(lb, blockSize)

      size_x = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
      size_y = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
      size_z = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1

      nxbBlock = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1
      nybBlock = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1
      nzbBlock = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1

      allocate(x(size_x))
      allocate(y(size_y))
      allocate(z(size_z))

      call Grid_getCellCoords(IAXIS, lb, CENTER, .true., x, size_x)
      call Grid_getCellCoords(JAXIS, lb, CENTER, .true., y, size_y)
      call Grid_getCellCoords(KAXIS, lb, CENTER, .true., z, size_z)

      ! Loop through cells, neglecting one layer of guard cells
      do ii = blkLimitsGC(LOW,IAXIS)+1, blkLimitsGC(HIGH,IAXIS)-1
         do jj = blkLimitsGC(LOW,JAXIS)+1, blkLimitsGC(HIGH,JAXIS)-1
            do kk = blkLimitsGC(LOW,KAXIS)+1, blkLimitsGC(HIGH,KAXIS)-1

               paccx = 0.0
               paccy = 0.0
               paccz = 0.0

               ! Loop over every sink particle
               do p = 1, localnpf

                  mass = particles_global(MASS_PART_PROP,p)

                  dx_inside = x(ii) - particles_global(POSX_PART_PROP,p)
                  dy_inside = y(jj) - particles_global(POSY_PART_PROP,p)
                  dz_inside = z(kk) - particles_global(POSZ_PART_PROP,p)

                  ! account for periodic BCs
                  do nx = -nrep_pbc, nrep_pbc
                     dx = dx_inside + nx*Lx
                     x2 = dx*dx
                     do ny = -nrep_pbc, nrep_pbc
                        dy = dy_inside + ny*Ly
                        y2 = dy*dy
                        do nz = -nrep_pbc, nrep_pbc
                           dz = dz_inside + nz*Lz
                           z2 = dz*dz

                           radius = sqrt(x2 + y2 + z2)

                           if (radius .lt. maxradius_pbc) then
                              if (radius .lt. softening_radius_comoving) then
                                 if (softeningtype.eq.1) then ! spline softening (see e.g., Price & Monaghan 2007)
                                    q = radius*hinv
                                    if ((q.gt.1.e-5).and.(q.lt.1.0)) &
                                       & kernelvalue = h2inv*(4.0/3.0*q-1.2*q**3+0.5*q**4)/radius
                                    if ((q.ge.1.0)  .and.(q.lt.2.0)) &
                                       & kernelvalue = h2inv * &
                                       & (8.0/3.0*q-3.0*q**2+1.2*q**3-1.0/6.0*q**4-1.0/(15.0*q**2))/radius
                                    paccx = paccx + dx*kernelvalue*mass
                                    paccy = paccy + dy*kernelvalue*mass
                                    paccz = paccz + dz*kernelvalue*mass
                                 endif
                                 if (softeningtype.eq.2) then ! linear kernel inside smoothing radius
                                    paccx = paccx + dx*slope*mass
                                    paccy = paccy + dy*slope*mass
                                    paccz = paccz + dz*slope*mass
                                 endif
                              else
                                 r3 = 1.0/radius**3
                                 paccx = paccx + dx*r3*mass
                                 paccy = paccy + dy*r3*mass
                                 paccz = paccz + dz*r3*mass
                              endif
                           endif ! within maxradius_pbc

                        enddo   ! nz
                     enddo   ! ny
                  enddo   ! nx

               end do   ! sinks

               ! x-acceleration:
               solnData(SGAX_VAR,ii,jj,kk) = paccx*prefactor

               ! y-acceleration:
               solnData(SGAY_VAR,ii,jj,kk) = paccy*prefactor

               ! z-acceleration:
               solnData(SGAZ_VAR,ii,jj,kk) = paccz*prefactor

            enddo   ! cells
         enddo   ! cells
      enddo   ! cells

      deallocate(x)
      deallocate(y)
      deallocate(z)

      call Grid_releaseBlkPtr(lb,solnData)

   end do   ! blocks

 end if   ! periodic BCs


 if (grav_boundary_type .eq. "isolated") then

   ! Loop through blocks
   do lb = 1, blockCount

      call Grid_getBlkPtr(blockList(lb),solnData)

      call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)
      call Grid_getBlkPhysicalSize(blockList(lb), blockSize)

      size_x = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
      size_y = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
      size_z = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1

      nxbBlock = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1
      nybBlock = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1
      nzbBlock = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1

      allocate(x(size_x))
      allocate(y(size_y))
      allocate(z(size_z))

      call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, .true., x, size_x)
      call Grid_getCellCoords(JAXIS, blockList(lb), CENTER, .true., y, size_y)
      call Grid_getCellCoords(KAXIS, blockList(lb), CENTER, .true., z, size_z)

      ! Loop through cells, neglecting one layer of guard cells
      do ii = blkLimitsGC(LOW,IAXIS)+1, blkLimitsGC(HIGH,IAXIS)-1
         do jj = blkLimitsGC(LOW,JAXIS)+1, blkLimitsGC(HIGH,JAXIS)-1
            do kk = blkLimitsGC(LOW,KAXIS)+1, blkLimitsGC(HIGH,KAXIS)-1

               paccx = 0.0
               paccy = 0.0
               paccz = 0.0

               ! Loop over global sink particles
               do p = 1, localnpf

                  mass = particles_global(MASS_PART_PROP,p)

                  dx = x(ii) - particles_global(POSX_PART_PROP,p)
                  dy = y(jj) - particles_global(POSY_PART_PROP,p)
                  dz = z(kk) - particles_global(POSZ_PART_PROP,p)

                  x2 = dx*dx
                  y2 = dy*dy
                  z2 = dz*dz

                  radius = sqrt(x2 + y2 + z2)

                  if (radius .lt. softening_radius_comoving) then
                     if (softeningtype.eq.1) then ! spline softening (see e.g., Price & Monaghan 2007)
                        q = radius*hinv
                        if ((q.gt.1.e-5).and.(q.lt.1.0)) &
                           & kernelvalue = h2inv*(4.0/3.0*q-1.2*q**3+0.5*q**4)/radius
                        if ((q.ge.1.0)  .and.(q.lt.2.0)) & 
                           & kernelvalue = h2inv * &
                           & (8.0/3.0*q-3.0*q**2+1.2*q**3-1.0/6.0*q**4-1.0/(15.0*q**2))/radius
                        paccx = paccx + dx*kernelvalue*mass
                        paccy = paccy + dy*kernelvalue*mass
                        paccz = paccz + dz*kernelvalue*mass
                     endif
                     if (softeningtype.eq.2) then ! linear kernel inside smoothing radius
                        paccx = paccx + dx*slope*mass
                        paccy = paccy + dy*slope*mass
                        paccz = paccz + dz*slope*mass
                     endif
                  else
                     r3 = 1.0/radius**3
                     paccx = paccx + dx*r3*mass
                     paccy = paccy + dy*r3*mass
                     paccz = paccz + dz*r3*mass
                  endif

               end do   ! sinks

               ! x-acceleration:
               solnData(SGAX_VAR,ii,jj,kk) = paccx*prefactor

               ! y-acceleration:
               solnData(SGAY_VAR,ii,jj,kk) = paccy*prefactor

               ! z-acceleration:
               solnData(SGAZ_VAR,ii,jj,kk) = paccz*prefactor

            enddo   ! cells
         enddo   ! cells
      enddo   ! cells

      deallocate(x)
      deallocate(y)
      deallocate(z)

      call Grid_releaseBlkPtr(blockList(lb),solnData)

   end do   ! blocks

 end if   ! isolated BCs

 return

end subroutine Particles_sinkAccelSinksOnGas
