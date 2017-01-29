!!****if* source/physics/Diffuse/DiffuseMain/Multigrid/diff_advanceTherm
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
!!      Supported boundary conditions are ??? and
!!      periodic (1).  The same boundary conditions are currently applied
!!      in all directions.
!!
!! ARGUMENTS
!!
!!   blockCount   : The number of blocks in the list
!!   blockList(:) : The list of blocks on which the solution must be updated
!!   dt           : The time step
!!   pass           : Ignored in unsplit solver.
!!                    pass=1 order of directional sweep X-Y-Z, 
!!                    pass=2 order of directional sweep Z-Y-X.
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
!!     WTMP_VAR:  contains source term that was passed to Grid_advanceDiffusion
!!     COND_VAR:  contains conductivity that was passed to Grid_advanceDiffusion
!!     DTMP_VAR:  contains temperature increment for the current simulation time.
!!  For the Multigrid implementation:
!!     ISLS_VAR (residual)
!!     ICOR_VAR (correction)
!!
!! NOTES
!!
!!  The interface of this subroutine must be explicitly know to code that
!!  calls it.  The simplest way to make it so is to have something like
!!     use diff_interface,ONLY: diff_advanceTherm
!!  in the calling routine.
!!***

!!REORDER(4): solnVec

subroutine diff_advanceTherm(blockCount,blockList,dt,pass)


  use Diffuse_data, ONLY : useDiffuse, diff_meshMe, diff_meshComm
  use diff_saData, ONLY : diff_boundary, &
       updateDiffuse, &
       diff_scaleFactThermSaTempDiff, diff_scaleFactThermSaTime
  use Eos_interface, ONLY : Eos_wrapped
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Conductivity_interface, ONLY : Conductivity
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
      Grid_advanceDiffusion, Grid_getBlkIndexLimits, Grid_fillGuardCells, &
      Grid_getDeltas

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Multigrid.h"
#include "Flash_mpi.h"

  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  real,intent(in) :: dt
  integer, OPTIONAL, intent(IN):: pass


  real, POINTER, DIMENSION(:,:,:,:) :: solnVec

  integer       :: ierr, i,j,k,n

  integer       :: lb
  integer       :: bcTypes(6)
  real          :: bcValues(2,6) = 0.
  integer       :: oldTemp
  real          :: chi
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits
  integer :: ncells
  real :: cond_zone, cond_avg, diff_coeff, xdens, xtemp
  real :: massfrac(NSPECIES)

  real, dimension(MDIM)          :: del
  integer                        :: numCells 
  real                           :: Cond_L, Cond_R, theta = 0.5
  logical                        :: solnIsDelta = .TRUE.

!=========================================================================
  if(.not.useDiffuse) return
  if(.not.updateDiffuse) return ! TO HANDLE DIFFUSION TYPE ?


  call Timers_start("diffusive advance Barrier")
  call MPI_Barrier (diff_meshComm, ierr)
  call Timers_stop("diffusive advance Barrier")

  call Timers_start("diffusive advance")

  bcTypes = MG_BND_PERIODIC ! For now 1-> Periodic, 2-> Dirichlet, see Multigrid.h
  bcValues = 0.

  call Grid_fillGuardCells(CENTER,ALLDIR)

  do lb = 1, blockCount

     call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)

     call Grid_getDeltas(blocklist(lb), del)

     call Grid_getBlkPtr(blocklist(lb), solnVec)

     do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
        do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
           do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)

                  xtemp = solnVec(TEMP_VAR,i,j,k)
                  xdens = solnVec(DENS_VAR,i,j,k)

                  ! load the mass fractions
                  do n = 1, NSPECIES 
                     massfrac(n) = solnVec(SPECIES_BEGIN-1+n,i,j,k)
                  enddo

                  call Conductivity(xtemp, xdens, massfrac, cond_zone, diff_coeff)

                  solnVec(COND_VAR,i,j,k) = diff_coeff
                  chi                     = diff_coeff

                  if (solnIsDelta == .TRUE.) then
                      solnVec(DTMP_VAR,i,j,k) = 0.0
                  else 
                      solnVec(DTMP_VAR,i,j,k) =  solnVec(TEMP_VAR,i,j,k)
                  endif
           end do
        end do
      end do

      

      do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
         do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
            do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

               Cond_R = chi !0.5*(solnVec(COND_VAR,i+1,j,k) + solnVec(COND_VAR,i,j,k))
               Cond_L = chi !0.5*(solnVec(COND_VAR,i-1,j,k) + solnVec(COND_VAR,i,j,k))

               if (solnIsDelta) then
                   SolnVec(WTMP_VAR,i,j,k) = (dt/del(1)**2) * &
                                             (cond_R*solnVec(TEMP_VAR,i+1,j,k) - (cond_R+cond_L)*solnVec(TEMP_VAR,i,j,k) + cond_L*solnVec(TEMP_VAR,i-1,j,k))
               else 
                   SolnVec(WTMP_VAR,i,j,k) = SolnVec(TEMP_VAR,i,j,k) + (1.0-theta)*(dt/del(1)**2) * &
                                            (cond_R*solnVec(TEMP_VAR,i+1,j,k) - (cond_R+cond_L)*solnVec(TEMP_VAR,i,j,k) + cond_L*solnVec(TEMP_VAR,i-1,j,k))
               endif


               if (NDIM .ge. 2) then
  
                   Cond_R = 0.5*(solnVec(COND_VAR,i,j+1,k) + solnVec(COND_VAR,i,j,k))
                   Cond_L = 0.5*(solnVec(COND_VAR,i,j-1,k) + solnVec(COND_VAR,i,j,k))

                   if (solnIsDelta) then
                       SolnVec(WTMP_VAR,i,j,k) = SolnVec(WTMP_VAR,i,j,k) + (dt/del(2)**2) * &
                                                (cond_R*solnVec(TEMP_VAR,i,j+1,k)-(cond_R+cond_L)*solnVec(TEMP_VAR,i,j,k)+cond_L*solnVec(TEMP_VAR,i,j-1,k))
                   else
                       SolnVec(WTMP_VAR,i,j,k) = SolnVec(WTMP_VAR,i,j,k) + (1.0-theta)*(dt/del(2)**2) * &
                                                (cond_R*solnVec(TEMP_VAR,i,j+1,k)-(cond_R+cond_L)*solnVec(TEMP_VAR,i,j,k)+cond_L*solnVec(TEMP_VAR,i,j-1,k))
                   endif
               endif

               if (NDIM .eq. 3) then

                   Cond_R = 0.5*(solnVec(COND_VAR,i,j,k+1) + solnVec(COND_VAR,i,j,k))
                   Cond_L = 0.5*(solnVec(COND_VAR,i,j,k-1) + solnVec(COND_VAR,i,j,k))

                   if (solnIsDelta) then
                       SolnVec(WTMP_VAR,i,j,k) = SolnVec(WTMP_VAR,i,j,k) + (dt/del(3)**2) * &
                                                (cond_R*solnVec(TEMP_VAR,i,j,k+1) - (cond_R+cond_L)*solnVec(TEMP_VAR,i,j,k) + cond_L*solnVec(TEMP_VAR,i,k-1,k))
                   else 
                       SolnVec(WTMP_VAR,i,j,k) =  SolnVec(WTMP_VAR,i,j,k) + (1.0-theta)*(dt/del(3)**2) * &
                                                 (cond_R*solnVec(TEMP_VAR,i,j,k+1) - (cond_R+cond_L)*solnVec(TEMP_VAR,i,j,k) + cond_L*solnVec(TEMP_VAR,i,k-1,k))
                   endif
              endif

           end do
        end do
     end do

     call Grid_releaseBlkPtr(blocklist(lb), solnVec)

  end do

   

   call Grid_advanceDiffusion(DTMP_VAR, WTMP_VAR, COND_VAR, -1, bcTypes, bcValues, &
                              dt* diff_scaleFactThermSaTime, chi, scaleFact=1.0,solnIsDelta=.TRUE.)
   

! WE HAVE NEW DTMP_VAR ... use it to update Solution

#ifdef USEBARS
  call MPI_Barrier (diff_meshComm, ierr)
#endif  
  call Timers_start("eos")
  do lb =1, blockCount
     call Grid_getBlkPtr(blocklist(lb), solnVec)

     if (solnIsDelta) then
         solnVec(TEMP_VAR,:,:,:) = solnVec(TEMP_VAR,:,:,:) + solnVec(DTMP_VAR,:,:,:)*diff_scaleFactThermSaTempDiff
     else
         solnVec(TEMP_VAR,:,:,:) = solnVec(DTMP_VAR,:,:,:)
     endif

     call Grid_releaseBlkPtr(blocklist(lb), solnVec)
  
     call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)

     call Eos_wrapped(MODE_DENS_TEMP, blkLimits, blockList(lb))
  enddo
      
  call Timers_stop("eos")

  call Timers_stop ("diffusive advance")
  
  return
end subroutine diff_advanceTherm
