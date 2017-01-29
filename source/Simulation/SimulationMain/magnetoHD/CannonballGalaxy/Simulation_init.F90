subroutine Simulation_init()

#include "Flash.h"

  use Simulation_data

#ifdef FLASH_USM_MHD
  use Hydro_data, ONLY : hy_prol_method !! HACK
#endif

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  use IO_interface, ONLY : IO_getScalar

  use Driver_interface, ONLY : Driver_abortFlash

  use Particles_data, ONLY : pt_posInitialized

  implicit none

#include "constants.h"
#include "Flash_mpi.h"
  
  ! Temporary variables.

  real, save :: kpc, Msun, keV, kms
  real       :: result1, result2, rmin
  integer :: ierr, i, numPEs, ii, jj, kk, totzones, halo

  real, external :: sim_getMagField
  logical :: partPosInitialized, updateRefine

  character(len=80) :: sim_prol_method_str
  character(len=1) :: dummy_char

!===============================================================================
  
  ! Initialization

  call mpi_comm_rank(MPI_COMM_WORLD, sim_meshMe, ierr)

  call PhysicalConstants_get("pi", sim_pi)

  sim_zSol = 0.02

  kpc = 3.0856775807E21
  Msun = 1.9889225E33
  k_B = 1.381e-16
  kms = 1.0e5
  N_a = 6.022e23
  mueinv = 0.875
  keV = 11604440.207109345
  m_e = 9.109e-28

  call RuntimeParameters_get("profile1", sim_profile1)
  call RuntimeParameters_get("profile2", sim_profile2)
  call RuntimeParameters_get("testSingleGalaxy", sim_testSingleGalaxy)
  call RuntimeParameters_get("testAtmosphere", sim_testAtmosphere)
  call RuntimeParameters_get("nsubzones", sim_subZones)
  call RuntimeParameters_get("xmax", sim_xMax)
  call RuntimeParameters_get("xmin", sim_xMin)
  call RuntimeParameters_get("ymax", sim_yMax)
  call RuntimeParameters_get("ymin", sim_ymin)
  call RuntimeParameters_get("zmax", sim_zMax)
  call RuntimeParameters_get("zmin", sim_zMin)
  call RuntimeParameters_get("d", sim_d)
  call RuntimeParameters_get("smlrho", sim_smlRho)
  call RuntimeParameters_get("smallx", sim_smallX)
  call RuntimeParameters_get("smallt", sim_smallT)
  call RuntimeParameters_get("smalle", sim_smalle)
  call RuntimeParameters_get("smallp", sim_smallP)
  call RuntimeParameters_get("RefinementDensityCutoff", &
       sim_refinementDensityCutoff)
  call RuntimeParameters_get("deRefiningRadius", sim_deRefiningRadius)
  call RuntimeParameters_get("RefiningRadius", sim_refiningRadius)
  call RuntimeParameters_get("restart", sim_restart)
  call RuntimeParameters_get("plasmaBeta", sim_plasmaBeta)
  call RuntimeParameters_get("lrefine_min", sim_lrefineMin)
  call RuntimeParameters_get("lrefine_max", sim_lrefineMax)
  call RuntimeParameters_get("useMeKaLCooling", sim_useMeKaLCooling)
  call RuntimeParameters_get("pressureNormalize", sim_pressureNormalize)
  call RuntimeParameters_get("Bmag", sim_Bmag)
  call RuntimeParameters_get("lMin", sim_lMin)
  call RuntimeParameters_get("lMax", sim_lMax)
  call RuntimeParameters_get("vInit", sim_vInit)
  call RuntimeParameters_get("ptdirn", sim_ptdirn)
#ifdef FLASH_USM_MHD
  call RuntimeParameters_get("prolMethod", sim_prol_method_str)
#endif

#ifdef MAGX_VAR
  call RuntimeParameters_get('killdivb',sim_killdivb)
#ifdef FLASH_USM_MHD
#include "UHD.h"
  call RuntimeParameters_get('forceHydroLimit',sim_forceHydroLimit)
  if(trim(sim_prol_method_str) == "injection_prol" .or. &
     trim(sim_prol_method_str) == "INJECTION_PROL" ) then
     hy_prol_method = INJECTION_PROL
  elseif (trim(sim_prol_method_str) == "balsara_prol" .or. &
          trim(sim_prol_method_str) == "BALSARA_PROL" ) then
     hy_prol_method = BALSARA_PROL
  else
     call Driver_abortFlash&
          ("[Simulation_init]: The prolongation method is of unknown type: " // &
           "Please choose one of 'injection_prol(2D & 3D)' or 'balsara_prol(3D)'.")
  endif
#endif
#else
  sim_killdivb = .false.
  sim_forceHydroLimit = .true.
#endif

  call mpi_comm_size(MPI_COMM_WORLD, numPEs, ierr)

  ! Compute derived quantities.

  call PhysicalConstants_get("Newton", sim_Newton)
  
  sim_d = sim_d*kpc
  sim_refiningRadius = sim_refiningRadius*kpc
  sim_deRefiningRadius = sim_deRefiningRadius*kpc
  sim_Bmag = sim_Bmag / sqrt(4.*PI)

  if (.not. sim_testSingleGalaxy) then

     if (sim_meshMe == 0) then
        
        open (9, file=sim_profile1)
        
        read (9, '(A,I7)') dummy_char, numPoints1
        
        print *, "numPoints1 = ", numPoints1
        
     endif
     
     call mpi_bcast(numPoints1, 1, FLASH_INTEGER, 0, MPI_COMM_WORLD, ierr)
     
     allocate(r1(numPoints1))
     allocate(dens1(numPoints1))
     allocate(pres1(numPoints1))
     allocate(gpot1(numPoints1))
     allocate(grav1(numPoints1))
     
     if (sim_meshMe == 0) then
        
        do i = 1, numPoints1
           
           read (9, *) r1(i), dens1(i), pres1(i), gpot1(i), grav1(i)
           
        enddo
        
        close(9)
        
     endif

     call mpi_bcast(r1, numPoints1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(dens1, numPoints1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(pres1, numPoints1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(gpot1, numPoints1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(grav1, numPoints1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)

  endif

  if (.not. sim_testAtmosphere) then

     if (sim_meshMe == 0) then
        
        open (10, file=sim_profile2)
        
        read (10, '(A,I7)') dummy_char, numPoints2
        
     endif

     call mpi_bcast(numPoints2, 1, FLASH_INTEGER, 0, MPI_COMM_WORLD, ierr)

     allocate(r2(numPoints2))
     allocate(dens2(numPoints2))
     allocate(pres2(numPoints2))
     allocate(gpot2(numPoints2))
     allocate(grav2(numPoints2))

     if (sim_meshMe == 0) then

        do i = 1, numPoints2
           
           read (10, '(5E14.7)') r2(i), dens2(i), pres2(i), gpot2(i), grav2(i)
        
        enddo

        close(10)

     endif

     call mpi_bcast(r2, numPoints2, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(dens2, numPoints2, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(pres2, numPoints2, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(gpot2, numPoints2, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(grav2, numPoints2, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)

  endif

  if (.not. sim_restart) then

     if (.not. sim_testAtmosphere) then

        if (sim_ptdirn == 1) then
           sim_xCtr = sim_d
           sim_yCtr = (sim_yMax + sim_yMin) / 2.
           sim_zCtr = (sim_zMax + sim_zMin) / 2.
        else if (sim_ptdirn == 2) then
           sim_xCtr = (sim_xMax + sim_xMin) / 2.
           sim_yCtr = sim_d
           sim_zCtr = (sim_zMax + sim_zMin) / 2.
        else if (sim_ptdirn == 3) then
           sim_xCtr = (sim_xMax + sim_xMin) / 2.
           sim_yCtr = (sim_yMax + sim_yMin) / 2.
           sim_zCtr = sim_d
        endif

        sim_rCtr = sim_d
        sim_vrCtr = sim_vInit
        sim_oarCtr = 0.
        
     endif

     nsubinv = 1./real(sim_subZones)
     nsubvolinv = nsubinv**3

     call mpi_barrier(MPI_COMM_WORLD, ierr)

     if (.not. sim_forceHydroLimit) then
        call sim_loadBfield(sim_meshMe)
     endif

  else

     if (.not. sim_testAtmosphere) then
       
        call IO_getScalar("subcluster r", sim_rCtr)
        call IO_getScalar("subcluster vr", sim_vrCtr)
        call IO_getScalar("subcluster ar", sim_oarCtr)

        if (sim_ptdirn == 1) then
           sim_xCtr = sim_rCtr
           sim_yCtr = (sim_yMax + sim_yMin) / 2.
           sim_zCtr = (sim_zMax + sim_zMin) / 2.
        else if (sim_ptdirn == 2) then
           sim_xCtr = (sim_xMax + sim_xMin) / 2.
           sim_yCtr = sim_rCtr
           sim_zCtr = (sim_zMax + sim_zMin) / 2.
        else if (sim_ptdirn == 3) then
           sim_xCtr = (sim_xMax + sim_xMin) / 2.
           sim_yCtr = (sim_yMax + sim_yMin) / 2.
           sim_zCtr = sim_rCtr
        endif

     endif

  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  
  return

end subroutine Simulation_init

subroutine sim_loadBfield()

  use Simulation_data

  implicit none

#include "Flash_mpi.h"
#include "Flash.h"

  integer :: ii, jj, kk, totzones, ierr

  lbox = NXB*(2**(sim_lrefineMin-1))
  
  sim_nBzones = lbox+NGUARD

  allocate(sim_Bxcoord(sim_nBzones))
  allocate(sim_Bycoord(sim_nBzones))
  allocate(sim_Bzcoord(sim_nBzones))

  sim_Bxmin = sim_xMin
  sim_Bymin = sim_yMin
  sim_Bzmin = sim_zMin
  sim_Bdx = (sim_xMax-sim_xMin)/lbox
  sim_Bdy = (sim_yMax-sim_yMin)/lbox
  sim_Bdz = (sim_zMax-sim_zMin)/lbox

  do ii = 1, sim_nBzones
     sim_Bxcoord(ii) = (real(ii)-1.0)*sim_Bdx + sim_Bxmin
  enddo

  do jj = 1, sim_nBzones
     sim_Bycoord(jj) = (real(jj)-1.0)*sim_Bdy + sim_Bymin
  enddo
  
  do kk = 1, sim_nBzones
     sim_Bzcoord(kk) = (real(kk)-1.0)*sim_Bdz + sim_Bzmin
  enddo

  allocate(sim_Bx(sim_nBzones,sim_nBzones,sim_nBzones))
  allocate(sim_By(sim_nBzones,sim_nBzones,sim_nBzones))
  allocate(sim_Bz(sim_nBzones,sim_nBzones,sim_nBzones))

  if (sim_meshMe == 0) call bsetup(sim_meshMe)
    
  totzones = sim_nBzones*sim_nBzones*sim_nBzones
  
  call mpi_bcast(sim_Bx, totzones, &
       FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(sim_By, totzones, &
       FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(sim_Bz, totzones, &
       FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  
  sim_BLx = sim_Bdx * sim_nBzones
  sim_BLy = sim_Bdy * sim_nBzones
  sim_BLz = sim_Bdz * sim_nBzones

  if (sim_meshMe == 0) print *, "Got the fields"

  return

end subroutine sim_loadBfield
