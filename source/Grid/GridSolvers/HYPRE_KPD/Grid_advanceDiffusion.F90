!!****if* source/Grid/GridSolvers/HYPRE_KPD/Grid_advanceDiffusion
!!
!! NAME
!!  Grid_advanceDiffusion
!!
!! SYNOPSIS
!!
!!  call Grid_advanceDiffusion (integer, intent(IN) :: iVar,
!!                              integer(IN)         :: iSrc,
!!                              integer, intent(IN) :: iFactorB,
!!                              integer, intent(IN) :: iFactorA,
!!                              integer, intent(IN) :: bcTypes(6),
!!                              real,    intent(IN) :: bcValues(2,6),
!!                              real,    intent(IN) :: dt,
!!                              real,    intent(IN) :: chi,
!!                              real,    intent(IN) :: scaleFact,
!!                              real,    intent(IN) :: theta,
!!                              integer, intent(IN), OPTIONAL :: iFactorC,
!!                              integer, intent(IN), OPTIONAL :: iFactorD
!!                              integer, OPTIONAL, intent(IN) :: pass,
!!                              logical(IN) :: solnIsDelta)
!!  DESCRIPTION 
!!
!!      This routine advances the diffusion operator of the form,
!!      A*(df/dt) + C*f = div(B*grad(f)) + D
!!      f -> Variable to be diffused.
!!      C,D are optional factors.
!!  
!!      Presently it is used to do conduction and multigroup diffusion.
!!  
!! ARGUMENTS
!!   iVar           : Variable on which the diffusion operatorion is performed (e.g TEMP_VAR)
!!   iFactorA       :| Are factors in the equation with spatial variation.
!!   iFactorB       :| Factor C,D are optional and are generally used
!!   iFactorC       :| to represent emission/absorption in MGD.
!!   iFactorD       :| iFactorA is needed only for conduction.
!!   bcTypes        : Presently OUTFLOW, VACUUM is supported, DIRICHLET is untested.
!!   bcValues       : Values of iVar,iFactorB on boundary (DIRICHLET).                        
!!   dt             : The time step.
!!   scaleFact      : Factor by which the end solution is scaled (not used).
!!   chi            : useful for constant diffusion problems (not used).
!!   theta          : varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank Nicholson
!!   pass           : Ignored in unsplit solver.
!!                    pass=1 order of directional sweep X-Y-Z, 
!!                    pass=2 order of directional sweep Z-Y-X.
!!   iSrc           : Ignored.
!!   solnIsDelta    : Is the solution only a delta that the caller has to apply to the
!!                    temperature, rather than temperature itself (ignored).
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!! Supports: 1D, 2D PARAMESH (with local refinement)
!!           1D, 2D, 3D in UG (3D untested).
!!           1D Spherical in PARAMESH/UG
!!           1D, 2D Cylindrical in PARAMESH/UG (R-Z)
!!
!!           Uses HYPRE library to solve AX = B
!!  
!!
!!***

subroutine Grid_advanceDiffusion (iVar, iSrc, iFactorB, iFactorA, bcTypes, bcValues, &
     dt, chi, scaleFact,theta, solnIsDelta,iFactorC, iFactorD, pass)
  
  use Grid_data,        ONLY : gr_meshMe, gr_meshcomm
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_interface,     ONLY : gr_hypreCreateMatrix, gr_hypreComputeB
  use Grid_interface,   ONLY : Grid_fillGuardCells, Grid_getListOfBlocks, &
                               Grid_getBlkPtr, Grid_releaseBlkPtr,        &
                               Grid_getBlkIndexLimits, Grid_getBlkData,   &
                               Grid_getBlkRefineLevel

  
  use gr_hypreData,   ONLY   : gr_hypreLower, gr_hypreUpper, &
                               gr_hypreMatA, gr_hypreVecB, gr_hypreRefineMIN

  implicit none
  
#include "Flash.h"
#include "constants.h"
  
  integer, intent(IN) :: iVar
  integer, intent(IN) :: iSrc
  integer, intent(IN) :: iFactorB
  integer, intent(IN) :: iFactorA
  real, intent(IN)    :: dt 
  real, intent(IN)    :: chi
  real, intent(IN)    :: scaleFact
  real, intent(IN)    :: theta
  logical, intent(IN) :: solnIsDelta
  integer, dimension(6),  intent(IN) :: bcTypes
  real   , dimension(2,6),intent(IN) :: bcValues
  integer, intent(IN), OPTIONAL :: pass
  integer, intent(IN), OPTIONAL :: iFactorC
  integer, intent(IN), OPTIONAL :: iFactorD   
  
  integer :: blockCount
  integer :: blockList(MAXBLOCKS)

  
!!! BEGIN TEMP

  integer :: mylevel, mypart, var, ii,i,j,k, ierr
  real, allocatable :: BoxVal(:)
  real, allocatable :: RHSVal(:)
  integer :: datasize(MDIM)
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer :: blockID, lb
  real, allocatable :: cellVolumes(:,:,:)
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits 
  
  character(len=32) :: matfile

!!! END TEMP
  
  call Timers_start("Grid_advanceDiffusion")   
  
!!$  call Timers_start("Diffusion Barrier")
!!$  call MPI_Barrier (gr_meshcomm, ierr)
!!$  call Timers_stop ("Diffusion Barrier")
  
  

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  
  !!-----------------------------------------------------------------------
  !!     1.  Do we need to reset HYPRE grid ?, has the underlying AMR
  !!         mesh been modified ?, is this the first call to HYPRE ?
  !!         Only in AMR is this routine meaningful.
  !!-----------------------------------------------------------------------
  call gr_hypreGridStatus (blockCount, blockList)  
  
  !!-----------------------------------------------------------------------
  !!     2. Exchange iFactorB across processors. Needs to be done only in
  !!        PARAMESH / AMR, if UG the function will return without any action.
  !!-----------------------------------------------------------------------
  call gr_hypreExchangeFacB (iFactorB, blockCount, blockList)  
  
  !!-----------------------------------------------------------------------
  !!     3. Set initial guess for solver typically 
  !!        iVar at previous time step is used.
  !!-----------------------------------------------------------------------
  call gr_hypreSetIniGuess (iVar, blockCount, blockList)  
    
  
  !!-----------------------------------------------------------------------
  !!     4. Compute the actual A matrix in AX=B
  !!-----------------------------------------------------------------------
  call gr_hypreCreateMatrix(iVar, iFactorB, iFactorA, bcTypes, bcValues, dt, theta,  &
       blockCount, blockList, .TRUE., iFactorC)     

  
  call gr_hypreComputeB (blockCount, blockList, iVar, iFactorA, iFactorB, dt, theta, &
       bcTypes, bcValues, iFactorD)  
  
  
  mypart = 0  !! part iterator 
  var    = 0  !! var iterator.
  
  do lb = 1, blockCount 
     
     blockID = blockList(lb)
     
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)    
     call Grid_getBlkPtr(blockID, solnVec)
     call Grid_getBlkRefineLevel(blockID,mylevel)

     mypart = mylevel - gr_hypreRefineMIN
     
     datasize(1:MDIM) = blkLimits(HIGH,1:MDIM)-blkLimits(LOW,1:MDIM)+1         
     
     allocate(cellVolumes(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))
     
     call Grid_getBlkData(blockID, CELL_VOLUME, 0, EXTERIOR,          &
          blkLimits(LOW,:), cellVolumes,datasize)      
     
     allocate(BoxVal(product(dataSize(1:NDIM))))
     allocate(RHSVal(product(dataSize(1:NDIM))))     
     
     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
              
              ii = (k - blkLimits(LOW,KAXIS)  + 1)                             +  &
                   (j - blkLimits(LOW,JAXIS))*dataSize(KAXIS)                  +  &
                   (i - blkLimits(LOW,IAXIS))*dataSize(KAXIS)*dataSize(JAXIS)  
              
              BoxVal(ii) = solnVec(iFactorA,i,j,k)*cellVolumes(i,j,k)
              
              RHSVal(ii) = BoxVal(ii)*solnVec(iVar,i,j,k)
              
              if (present(iFactorC)) then
                 BoxVal(ii) = BoxVal(ii) + dt*solnVec(iFactorC,i,j,k)*cellVolumes(i,j,k)         
              end if
              
              if (present(iFactorD)) then 
                 RHSVal(ii) =  RHSVal(ii) + SolnVec(iFactorD,i,j,k)*cellVolumes(i,j,k)*dt
              end if
                            
           end do
        end do
     end do
     
     
     call HYPRE_SStructVectorAddToBoxValu(gr_hypreVecB, mypart, gr_hypreLower(lb, 1:NDIM), &
          gr_hypreUpper(lb,1:NDIM), var, RHSVal(:), ierr)
     
     
     call HYPRE_SStructMatrixAddToBoxValu(gr_hypreMatA, mypart, gr_hypreLower(lb, 1:NDIM), & 
          gr_hypreUpper(lb,1:NDIM), var, 1, (/0/), BoxVal(:), ierr)     
     
     
     call Grid_releaseBlkPtr(blockID, solnVec)
     
     deallocate (BoxVal)
     deallocate (RHSVal)
     deallocate(cellVolumes)
     
  end do  
  
  
!!$  matfile = 'ex12f.out'
!!$  matfile(10:10) = char(0)
!!$  call HYPRE_SStructMatrixPrint(matfile, gr_hypreMatA, 0, ierr)
!!$  pause
  
  
  !!-----------------------------------------------------------------------
  !!     7. Solve AX = B
  !!-----------------------------------------------------------------------
  call gr_hypreSolve ()
  
  !!-----------------------------------------------------------------------
  !!     8. Update unk variable using X.
  !!-----------------------------------------------------------------------
  call gr_hypreUpdateSoln (iVar, blockCount, blockList)  
  
  
  call Timers_stop("Grid_advanceDiffusion") 
  
  return
  
end subroutine Grid_advanceDiffusion
