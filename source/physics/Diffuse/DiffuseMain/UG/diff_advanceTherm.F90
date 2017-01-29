!!****if* source/physics/Diffuse/DiffuseMain/UG/diff_advanceTherm
!!
!!  NAME 
!!
!!  diff_advanceTherm
!!
!!  SYNOPSIS
!!
!!  call diff_advanceTherm(integer(IN)                  :: blockCount,
!!                         integer(IN)                  :: blockList(blockCount),
!!                         real(IN)                     :: dt,
!!                         integer, OPTIONAL, intent(IN):: pass)
!!
!!  DESCRIPTION 
!!      This routine advances the heat diffusion equation (i.e., heat conduction).
!!      An implicit scheme is used. 
!!
!!      Supported boundary conditions are: 
!!                PERIODIC, OUTFLOW (tested).
!!                DIRICHLET (untested).
!!
!!
!! ARGUMENTS
!!
!!  blockCount   - The number of blocks in the list
!!  blockList(:) - The list of blocks on which the solution must be updated
!!   dt           : The time step
!!   pass         : pass=1 directional order of solution sweep X-Y-Z, 
!!                  pass=2 directional order of solution sweep Z-Y-X.
!!  dt           - The time step
!!
!! SIDE EFFECTS
!!
!!  Updates certain variables in permanent UNK storage to contain the
!!  updated temperature and some auxiliaries.  Invokes a solver (of the Heat
!!  diffusion equation). On return,
!!     TEMP_VAR:  contains updated temperature for the current simulation time.
!!     EINT_VAR, ENER_VAR, PRES_VAR, etc.: updated accordingly by EOS call
!!
!!  May modify certain variables used for intermediate results by the solvers
!!  invoked. The list of variables depends on the Diffuse implementation.
!!  The following information is subject to change.
!!     COND_VAR:  contains conductivity that was passed to Grid_advanceDiffusion
!!  
!!  NOTES:
!!  
!!  The current implementation of Radiation diffusion is a Gray Approximation.
!!  Types: P1, Simple Diffusion, Flux limited Diffusion, P1/3 and Variable Eddington Factor (VEF)
!!  We have implemented Simple Diffusion.
!!
!!
!! NOTES
!!
!!  The interface of this subroutine must be explicitly know to code that
!!  calls it.  The simplest way to make it so is to have something like
!!  use diff_interface,ONLY: diff_advanceTherm
!!  in the calling routine.
!!***

!!REORDER(4): solnVec

subroutine diff_advanceTherm(blockCount,blockList,dt,pass)


  use Diffuse_data, ONLY : useDiffuse, diff_meshMe, diff_meshcomm,&
       diff_useRadDiffusion, diff_useEleCond, diff_useERad,diff_asol, &
       diff_useRadFlxLimiter, diff_mele, diff_boltz, &
       diff_singleSpeciesA, diff_singleSpeciesZ, diff_avo, diff_mele, &
       diff_eleFlMode, diff_eleFlCoef
  
  use diff_saData, ONLY : diff_boundary, &
       updateDiffuse, &
       diff_scaleFactThermSaTempDiff, diff_scaleFactThermSaTime, &
       diff_eleDomainBC, diff_radDomainBC
  use Eos_interface, ONLY : Eos_wrapped, Eos, Eos_getAbarZbar
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Conductivity_interface, ONLY : Conductivity
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
      Grid_advanceDiffusion, Grid_getBlkIndexLimits, Grid_fillGuardCells, &
      Grid_getDeltas, GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
      GRID_PDE_BND_DIRICHLET
  use Diffuse_interface, ONLY: Diffuse_solveScalar, &
       Diffuse_fluxLimiter

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"
#include "Eos.h"

!#define USE_ERAD_VAR 0

  integer,intent(IN)                       :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  real,intent(in)                          :: dt
  integer, OPTIONAL, intent(IN)            :: pass

  real, POINTER, DIMENSION(:,:,:,:) :: solnVec

  integer       :: oldTemp
  integer       :: ierr, i,j,k,n
  integer       :: lb
  integer       :: bcTypes(6)

  real          :: outputScaleFact
  real          :: bcValues(2,6) = 0.
  integer       :: Temptodiffuse, cvToUse
  real          :: chi
  real          :: cond_zone, diff_coeff, xdens, xtemp
  real          :: massfrac(NSPECIES), Ye
  real          :: gradE(NDIM), R, D

  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits
  real, dimension(MDIM) :: del

  ! To compute CV
  real,    dimension(EOS_NUM)            :: eos_arr
  logical, dimension(EOS_VARS+1:EOS_NUM) :: mask
  integer                                :: mode, vecLen
  logical                                :: changedZone
  logical                                :: ldiff


!=========================================================================
  if(.not.useDiffuse) return
  if(.not.updateDiffuse) return

  changedZone = .FALSE.

  call Timers_start("diffusive advance Barrier")
  call MPI_Barrier (diff_meshcomm, ierr)
  call Timers_stop("diffusive advance Barrier")

  call Timers_start("diffusive advance")

  outputScaleFact=1.0

  call Grid_fillGuardCells(CENTER,ALLDIR)

#ifdef TRAD_VAR

  !Radiation Diffusion - Setup BC Type
  bcTypes(:) = diff_radDomainBC(:)
  bcValues = 0.

  ! Note: We have two options 1) Diffuse using ERAD or 2) Diffuse using TRAD
  ! Since we store eRAD  (= ERAD / DENS), when we do the diffusion operation using ERAD
  ! we multiply eRAD (stored variable) by DENS and then later dived by DENS.

  where (bcTypes == PERIODIC)
     bcTypes = GRID_PDE_BND_PERIODIC
  elsewhere (bcTypes == DIRICHLET)
     bcTypes = GRID_PDE_BND_DIRICHLET
  elsewhere (bcTypes == OUTFLOW)
     bcTypes = GRID_PDE_BND_NEUMANN
  end where

  if (minval(bcTypes) == -1) then
      call Driver_abortFlash("[diff_advanceTherm] Boundary type not supported.")
  end if

! Implicit heat eq advance is solved with the total oldTemp of TEMP_VAR
  oldTemp=TEMP_VAR

  ! Radiation diffusion equation
  ! dEr/dt = d/dx (c/3k*dEr/dx) => DFCF_VAR = -1 & COND_VAR = c/3k
  ! where Er = ERAD_VAR * rho
  ! Er = aT^4 => dTr/dt = (1/4aTr^3)*d/dx(4aTr^3*c/3k*dTr/dx) => DFCF_VAR=1/4aT^3 & COND_VAR = 4aT^3*c/3k

  if (diff_useRadDiffusion) then !! PowerLaw Radiation Diffusion.
 
     do lb = 1, blockCount
        call Grid_getDeltas(blocklist(lb), del)
        call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
        call Grid_getBlkPtr(blocklist(lb), solnVec)
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 xtemp = solnVec(TRAD_VAR,i,j,k)
                 xdens = solnVec(DENS_VAR,i,j,k)
                 ! load the mass fractions
                 do n = 1, NSPECIES
                    massfrac(n) = solnVec(SPECIES_BEGIN-1+n,i,j,k)
                 enddo
                 call Conductivity(xtemp, xdens, massfrac, cond_zone, diff_coeff, 3)      ! 3: Rad Diff, 2: Elec Diff
                 if (diff_useERad) then
                     solnVec(ERAD_VAR,i,j,k) = solnVec(ERAD_VAR,i,j,k)*solnVec(DENS_VAR,i,j,k)

                     if (diff_useRadFlxLimiter) then 
                         gradE(1) = abs (solnVec(ERAD_VAR,i+1,j,k) - solnVec(ERAD_VAR,i-1,j,k))/del(1)
                         R        = gradE(1) / (cond_zone*solnVec(ERAD_VAR,i,j,k))
                         D        = 3.0 / (3.0 + R)
                        solnVec(COND_VAR,i,j,k) = D*diff_coeff   
                     else
                        solnVec(COND_VAR,i,j,k) = diff_coeff 
                     end if
                     solnVec(DFCF_VAR,i,j,k) = 1.0
                     
                 else
                     ! What if xtemp = 0.0 ???
                     solnVec(COND_VAR,i,j,k) = diff_coeff*(4*diff_asol*xtemp**3)
                     solnVec(DFCF_VAR,i,j,k) = 4*diff_asol*xtemp**3
                 endif
                 chi = diff_coeff

              enddo
           enddo
        enddo
        call Grid_releaseBlkPtr(blocklist(lb), solnVec)
     end do

     if (diff_useERad) then
 
         call Diffuse_solveScalar(ERAD_VAR, COND_VAR, DFCF_VAR,                   &
                                  bcTypes, bcValues, dt, outputScaleFact,         &
                                  chi, 0.5, pass, blockCount,blockList)               

         do lb = 1, blockCount
            call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
            call Grid_getBlkPtr(blocklist(lb), solnVec)
            do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
               do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
                  do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
                     solnVec (ERAD_VAR,i,j,k) = (solnVec(ERAD_VAR,i,j,k)/solnVec(DENS_VAR,i,j,k)) 
                  enddo
               enddo
           enddo
           ! What we SHOULD use here is something like call Eos_wrapped(MODE_DENS_EI_RAD,...), if that worked. - KW
           call Eos_wrapped(MODE_DENS_EI_GATHER, blkLimits, blockList(lb)) ! Update EOS
           call Grid_releaseBlkPtr(blocklist(lb), solnVec)
        end do


        
     else
        call Diffuse_solveScalar(TRAD_VAR, COND_VAR, DFCF_VAR,                &
                                 bcTypes, bcValues, dt, outputScaleFact,      &
                                 chi, 0.5, pass, blockCount,blockList)

        changedZone = .TRUE.
     endif

  endif
   
#endif



  ldiff = .false.
#ifdef TELE_VAR
  if (diff_useEleCond) then 
     Temptodiffuse = TELE_VAR 
     cvToUse       =  EOS_CVELE
     ldiff         = .true.

     !Electron Conduction - Setup BC Type
     bcTypes(:) = diff_eleDomainBC(:)
     bcValues = 0.
     where (bcTypes == PERIODIC)
        bcTypes = GRID_PDE_BND_PERIODIC
     elsewhere (bcTypes == DIRICHLET)
        bcTypes = GRID_PDE_BND_DIRICHLET
     elsewhere (bcTypes == OUTFLOW)
        bcTypes = GRID_PDE_BND_NEUMANN
     end where
  end if
#else
  if (useDiffuse) then
     Temptodiffuse=TEMP_VAR ! The default would follow Electron power law.
     cvToUse = EOS_CV
     ldiff = .true.

     !1T - setup
     bcTypes(:) = diff_eleDomainBC(:)
     bcValues = 0.
     where (bcTypes == PERIODIC)
        bcTypes = GRID_PDE_BND_PERIODIC
     elsewhere (bcTypes == DIRICHLET)
        bcTypes = GRID_PDE_BND_DIRICHLET
     elsewhere (bcTypes == OUTFLOW)
        bcTypes = GRID_PDE_BND_NEUMANN
     end where
  end if
#endif



  ! cv*rho*dT/dt = d/dx(k*dT/dx) => COND_VAR=k & DFCF_VAR=cv*rho
  if(ldiff) then
     do lb = 1, blockCount
        call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
        call Grid_getBlkPtr(blocklist(lb), solnVec)
        do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
           do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
              do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
                 xtemp = solnVec(Temptodiffuse,i,j,k)
                 xdens = solnVec(DENS_VAR,i,j,k)
                 ! load the mass fractions
                 do n = 1, NSPECIES
                    massfrac(n) = solnVec(SPECIES_BEGIN-1+n,i,j,k)
                 enddo
                 call Conductivity(xtemp, xdens, massfrac, cond_zone, diff_coeff, 2)
                 
                 ! Compute CV
                 vecLen = 1
                 mode = MODE_DENS_TEMP
                 eos_arr(Temptodiffuse) = xtemp
                 eos_arr(EOS_DENS) = xdens
                 mask = .false.
                 mask(cvToUse) = .true.
                 mask(EOS_DET) = .true.

                 call Eos(mode,vecLen,eos_arr,massfrac,mask)

                 solnVec(COND_VAR,i,j,k) = cond_zone
                 solnVec(DFCF_VAR,i,j,k) = xdens*eos_arr(cvToUse)
                 
                 ! Set abar and zbar:
#ifdef FLASH_MULTISPECIES
                 call Eos_getAbarZbar(solnVec(:,i,j,k),Ye=Ye,massFrac=massfrac)
#else
#ifdef YE_MSCALAR
                 Ye = solnVec(YE_MSCALAR,i,j,k)
#else
                 Ye = diff_singleSpeciesZ / diff_singleSpeciesA
#endif
#endif           
                 ! Set electron flux limiter:
                 solnVec(FLLM_VAR,i,j,k) = diff_eleFlCoef * &
                      sqrt(diff_boltz*xtemp/diff_mele) * &
                      diff_boltz*xtemp * &
                      (Ye * diff_avo * xdens)
              enddo
           enddo
        enddo
        
        call Grid_releaseBlkPtr(blocklist(lb), solnVec)
     end do
     
     call Diffuse_fluxLimiter(COND_VAR, Temptodiffuse, FLLM_VAR, &
          diff_eleFlMode, blockCount, blockList)
     
     call Diffuse_solveScalar(Temptodiffuse, COND_VAR, DFCF_VAR, bcTypes, &
          bcValues, dt, outputScaleFact, chi, 0.5, pass, blockCount,blockList)
     
     
     
     changedZone = .TRUE.
  endif
#ifdef USEBARS
  call MPI_Barrier (diff_meshcomm, ierr)
#endif 
  
  if (changedZone) then 
     call Timers_start("eos")
     do lb =1, blockCount
        call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)              
#ifdef TELE_VAR
        call Eos_wrapped(MODE_DENS_TEMP_GATHER, blkLimits, blockList(lb))
#else
        call Eos_wrapped(MODE_DENS_TEMP, blkLimits, blockList(lb))
#endif 
     enddo
      
     call Timers_stop("eos")
  endif
  
  call Timers_stop ("diffusive advance")
  
  return
end subroutine diff_advanceTherm
