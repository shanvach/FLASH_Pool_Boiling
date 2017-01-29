!!****if* source/Grid/GridSolvers/HYPRE_KPDa/paramesh/gr_hypreCreateMatrix_KPD
!!
!!  NAME 
!!
!!  gr_hypreCreateMatrix
!!
!!  SYNOPSIS
!!
!!  call gr_hypreCreateMatrix (integer, intent(IN)           :: iVar,
!!                             integer, intent(IN)           :: iFactorB,
!!                             integer, intent(IN)           :: iFactorA,
!!                             integer, intent(IN)           :: bcTypes(6),
!!                             real,    intent(IN)           :: bcValues(2,6),
!!                             real,    intent(IN)           :: dt,
!!                             real,    intent(IN)           :: alpha,
!!                             integer,           intent(IN) :: blockCount,
!!                             integer,dimension(blockCount),intent(IN) :: blockList,
!!                             logical, intent(IN)           :: JacobiMatrix
!!                             integer, intent(IN), OPTIONAL :: iFactorC)
!!
!!
!!  DESCRIPTION 
!!      This routine computes two matrices
!!          AX = B, where A is the matrix to be inverted
!!          B = MX, where M is a matrix whose product with iVar produces RHS B.
!!          if present iFactorC is added to the diagnol.
!!
!!      A*(df/dt) + C*f = div(B*grad(f)) + D
!!      f -> Variable to be diffused.
!!      C,D are optional factors.
!!
!!
!! ARGUMENTS
!! 
!!   iVar         : Variable on which the diffusion operatorion is performed (e.g TEMP_VAR)
!!   iFactorB     : Are factors in the equation with spatial variation. Factor C,D are optional 
!!   iFactorA     : and are generally used to represent emission/absorption in MGD.
!!   iFactorC     :
!!   bcTypes      : Presently OUTFLOW, VACUUM is supported, DIRICHLET is untested.
!!   bcValues     : Values of iVar,iFactorB on boundary (DIRICHLET, not used).                        
!!   dt           : The time step.
!!   alpha        : varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank Nicholson
!!   blockCount   : The number of blocks in the list.   
!!   blockList    : The list of blocks on which the solution must be updated.        
!!   JacobiMatrix : TRUE, computes A and FALSE computes M.
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!  
!!
!!***

!!REORDER(4): solnVec

!subroutine gr_hypreCreateMatrix(iVar, iFactorB, iFactorA, bcTypes, bcValues, dt, &
!     alpha, blockCount, blockList, JacobiMatrix, iFactorC, level, idenvar)

subroutine gr_hypreCreateMatrix_KPD(iVar, bcTypes, bcValues, dt, &
     alpha, blockCount, blockList, JacobiMatrix, level)
  
  use gr_hypreData,     ONLY : gr_hypreSolverType,gr_hypreLower, gr_hypreUpper, &
                               gr_hypreMatA, gr_hypreSetup, &
                               gr_hypreRefineMIN, gr_hypreRefineMAX, gr_hypreNeghLevels, &
                               gr_asol, gr_speedlt , gr_hypreSurrBlkSum
  use Grid_interface,   ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
                               Grid_getBlkIndexLimits, Grid_fillGuardCells, &
                               Grid_getBlkBC, Grid_getBlkCornerID,  &
                               Grid_getCellCoords, Grid_getFluxData, Grid_getBlkData, &
                               Grid_getLocalNumBlks
  use Timers_interface, ONLY : Timers_start, Timers_stop 
  use Grid_interface,   ONLY : GRID_PDE_BND_PERIODIC,  &
                               GRID_PDE_BND_NEUMANN,   &
                               GRID_PDE_BND_DIRICHLET
  use Grid_data,        ONLY : gr_geometry, gr_meshMe
  
  use tree, only : lrefine
  
  implicit none
 
#include "Flash.h" 
#include "constants.h"
#include "Flash_mpi.h"
#include "HYPREf.h"    
  
  !!-----------------------------------------------------------------------
  !!         ARGUMENTS
  !!-----------------------------------------------------------------------
  integer, intent(IN) :: iVar, level
  integer, intent(IN) :: bcTypes(6)
  real,    intent(IN) :: bcValues(2,6)
  real,    intent(IN) :: dt
  real,    intent(IN) :: alpha
  integer, intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  logical, intent(IN) :: JacobiMatrix
      
  !!-----------------------------------------------------------------------
  !!         LOCAL VARIABLES.
  !!-----------------------------------------------------------------------  
  integer :: ierr, pos(NDIM)
  real, dimension(MDIM)     :: del , delph, delmh
  real, dimension(2*MDIM, MDIM) :: negh_del
  real :: CDiv 
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits 
  integer :: datasize(MDIM), datasizeGC(MDIM)
  integer ::  mypart,mylevel
  integer ::  var
  integer :: blockID, lnblocks
  integer ::  nentries,  stencil_indices(7)  
  real    ::  values(7), graph_value(7) 
  integer :: i, j, k, ii, lb, lb_AMR
!!$  integer :: left, center, right
  real   :: condiph, condimh
  real   :: condjph, condjmh
  real   :: condkph, condkmh
  real   :: coeff, theta
  integer, dimension(2,MDIM):: faces 
  integer :: eachNegh,numNegh 

!  real, allocatable, dimension(:,:,:,:) :: xflux, yflux, zflux
  real :: dirichlet_multiplier
  integer :: dir
!  real, allocatable :: faceAreas  (:,:,:,:)  
  character(len=32) :: matfile
  integer :: numGraph, iter
!  real, allocatable :: BoxVal(:)
  real, dimension(28672) :: temp_BoxVal
  integer, parameter :: nFluxVars = 2**(NDIM-1)

  !- kpd - 
  real :: MdensXL, MdensXR, MdensYL, MdensYR, MdensZL, MdensZR
  real, pointer, dimension(:,:,:,:) :: facevarx,facevary,facevarz,solnData 
  integer :: iBoxCount

  call Timers_start("gr_hypreCreateMatrix")  
  
  if(JacobiMatrix) then
     theta = alpha
     dirichlet_multiplier = 1.0
  else
     theta = alpha - 1.0
     dirichlet_multiplier = 0.0     
  end if
  
  nentries = 2*NDIM + 1 
  do i = 1, nentries
     stencil_indices(i) = i-1
  enddo
  
  mypart = 0 !! part iterator 
  var    = 0  !! var iterator.  
  
  if (blockCount > 0) then
     call Grid_getBlkIndexLimits(blockList(1),blkLimits,blkLimitsGC)         
     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1          
!     allocate(xflux(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
!     allocate(yflux(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
!     allocate(zflux(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
  else
     datasize = 0.
  end if
  


  if(blockCount == 0) then
     call Timers_start("gr_hypreApplyBcToEdge")     
     call Timers_stop("gr_hypreApplyBcToEdge")     
  end if
  
  call Grid_getLocalNumBlks(lnblocks)

  !****************************************************************************
  !****************************************************************************
  iBoxCount = 0
  !do lb = 1, blockCount     
  !   blockID = blockList(lb)
  lb = 0
  do lb_AMR = 1, lnblocks
     if (lrefine(lb_AMR) == level) then

     blockID = lb_AMR
     lb = lb+1

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)    
     call Grid_getBlkPtr(blockID, solnVec)
     call Grid_getDeltas(blockID, del)
     call Grid_getBlkBC (blockID, faces)     
     call Grid_getBlkRefineLevel(blockID,mylevel)
    
     call Grid_getBlkPtr(blockID,facevarx ,FACEX)
     call Grid_getBlkPtr(blockID,facevary ,FACEY)
#if NDIM == 3
     call Grid_getBlkPtr(blockID,facevarz ,FACEZ)
#endif     
 
     datasize  (1:MDIM)= blkLimits  (HIGH,1:MDIM)-blkLimits  (LOW,1:MDIM)+1
     datasizeGC(1:MDIM)= blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1
     
!     allocate(BoxVal(nentries*product(datasize(1:NDIM))))                
     
     mypart = mylevel - gr_hypreRefineMIN                  
     
!     allocate(faceAreas(NDIM, blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),   &
!                              blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),   &
!                              blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))               
     
!     call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR, & 
!          blkLimitsGC(LOW,:), faceAreas(IAXIS,:,:,:), datasizeGC)         
     
#if NDIM >= 2
     
!     call Grid_getBlkData(blockID, CELL_FACEAREA, JLO_FACE, EXTERIOR, &
!          blkLimitsGC(LOW,:), faceAreas(JAXIS,:,:,:), datasizeGC)        
     
#if NDIM == 3
     
!     call Grid_getBlkData(blockID, CELL_FACEAREA, KLO_FACE, EXTERIOR, &
!          blkLimitsGC(LOW,:), faceAreas(KAXIS,:,:,:), datasizeGC)        

#endif     
#endif          


     if (mylevel > gr_hypreNeghLevels(lb,LEFT_EDGE,1+K2D,1+K3D)) then !!F/C
        negh_del(1,:) = del(:)*2.0
     else !! C/F
        negh_del(1,:) = del(:)/2.0
     end if
     
     if (mylevel > gr_hypreNeghLevels(lb,RIGHT_EDGE,1+K2D,1+K3D)) then !!F/C
        negh_del(2,:) = del(:)*2.0
     else !! C/F
        negh_del(2,:) = del(:)/2.0
     end if
     
     if (mylevel > gr_hypreNeghLevels(lb,1+K2D,LEFT_EDGE,1+K3D)) then !!F/C
        negh_del(3,:) = del(:)*2.0
     else !! C/F
        negh_del(3,:) = del(:)/2.0
     end if
     
     if (mylevel > gr_hypreNeghLevels(lb,1+K2D,RIGHT_EDGE,1+K3D)) then !!F/C
        negh_del(4,:) = del(:)*2.0
     else !! C/F
        negh_del(4,:) = del(:)/2.0
     end if
     
     if (mylevel > gr_hypreNeghLevels(lb,1+K2D,1+K2D,LEFT_EDGE)) then !!F/C
        negh_del(5,:) = del(:)*2.0
     else !! C/F
        negh_del(5,:) = del(:)/2.0
     end if
     
     if (mylevel > gr_hypreNeghLevels(lb,1+K2D,1+K2D,RIGHT_EDGE)) then !!F/C
        negh_del(6,:) = del(:)*2.0
     else !! C/F
        negh_del(6,:) = del(:)/2.0
     end if


!     call Grid_getFluxData(blockID,IAXIS,xflux,datasizeGC)
#if NDIM >= 2
!     call Grid_getFluxData(blockID,JAXIS,yflux,datasizeGC)
#if NDIM == 3
!     call Grid_getFluxData(blockID,KAXIS,zflux,datasizeGC)
#endif
#endif     
     
     iter = 1
     
     !BoxVal = 0.0
     temp_BoxVal = 0.0
     
     !print*,"Mdens is hardwired to 1 in createMatrix_KPD"

     do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)                                   
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)                    
              
              !!-------------------------------------------
              !!- kpd - For variable density implementation
              !!-------------------------------------------
              !MdensXL = facevarx(idenvar,i  ,j  ,1)
              !MdensXR = facevarx(idenvar,i+1,j  ,1)
              !MdensYL = facevary(idenvar,i  ,j  ,1)
              !MdensYR = facevary(idenvar,i  ,j+1,1)

              !!------------------------------------------------------------------------
              !!- kpd - Safeguard for guard cell filling and restriction/prolongation...
              !!------------------------------------------------------------------------
              !if (MdensXL==0. .OR. MdensXR==0. .OR. MdensYL==0. .OR. MdensYR==0. ) then
              !   print*,"================================================================================"
              !   print*,"ERROR hypreCreate: Density Equals Zero at a FACE",myPE,lb,i,j,MdensXR,MdensYR
              !   print*,"================================================================================"
              !elseif (MdensXL.lt.0.001 .OR. MdensXR.lt.0.001 .OR. MdensYL.lt.0.001 .OR. MdensYR.lt.0.001 ) then
              !   print*,"================================================================================"
              !   print*,"ERROR hypreCreate: Inverse Density Less than 0.001",myPE,lb,i,j,MdensXR,MdensYR
              !   print*,"================================================================================"
              !elseif (MdensXL.gt.1.0 .OR. MdensXR.gt.1.0 .OR. MdensYL.gt.1.0 .OR. MdensYR.gt.1.0 ) then
              !   print*,"================================================================================"
              !   print*,"ERROR hypreCreate: Inverse Density Greater than 1.0",myPE,lb,i,j,MdensXR,MdensYR
              !   print*,"================================================================================"
              !end if

              !!-----------------------------------------------------------------------
              !!         SWAP INDICES FOR FLASH-HYPRE COMPATIABILITY.
              !!-----------------------------------------------------------------------
#if NDIM == 1
              pos(1) = gr_hypreLower(lb,1) + (i-blkLimits(LOW,IAXIS))
#endif
#if NDIM == 2                                 
              pos(1) = gr_hypreLower(lb,1) + (j-blkLimits(LOW,JAXIS))
              pos(2) = gr_hypreLower(lb,2) + (i-blkLimits(LOW,IAXIS))
#endif
#if NDIM ==3
              pos(1) = gr_hypreLower(lb,1) + (k-blkLimits(LOW,KAXIS))
              pos(2) = gr_hypreLower(lb,2) + (j-blkLimits(LOW,JAXIS))
              pos(3) = gr_hypreLower(lb,3) + (i-blkLimits(LOW,IAXIS))
#endif         
              
              !!-----------------------------------------------------------------------
              !!         STENCILED VALUES TO BE FED INTO HYPRE MATRIX.
              !!-----------------------------------------------------------------------
              values = 0.0                           
              
              delmh = del
              delph = del
              
              condimh = 0.
              condiph = 0.
              
              numGraph = 0
              
              !! i-1,j,k
              if ((i /= blkLimits(LOW, IAXIS)) .or. (faces(1,IAXIS) == NOT_BOUNDARY)) then                                         

                 MdensXL = facevarx(MGW8_FACE_VAR,i  ,j  ,k)              
                 !MdensXL = facevarx(RH1F_FACE_VAR,i  ,j  ,k) +facevarx(RH2F_FACE_VAR,i  ,j  ,k)             
                 !print*,"Mdens is hardwired to 1 in createMatrix_KPD"
                 !MdensXL = 1.0 

                 if (MdensXL .lt. 0.99 .OR. MdensXL .gt. 1300.0) then 
                    print*,"ERROR in hypreCreateMatrix: Density is out of physical bounds!"
                 end if

                 if (i == blkLimits(LOW, IAXIS) .and. mylevel /= gr_hypreNeghLevels(lb,LEFT_EDGE, 1+K2D, 1+K3D)) then  !! F/C boundary.                                         
                 !   - kpd - Put in Driver Abort for fine/coarse interface in HYPRE

                 else                    

                    !! STENCILED RELATIONSHIP.                       
                    !*********************************************
                    !- kpd - For FLASH4 implementation for i-1,j,k
                    !*********************************************
                    temp_BoxVal(iter+1) = -1.d0 / (del(IAXIS)**2.0) * MdensXL 
      
                    !- kpd - Normalize the pressure field
                    if (lb .eq. 1 .AND. gr_meshMe .eq. 0 .AND.             &
                                        i .eq. blkLimits(LOW, IAXIS) .AND. &
                                        j .eq. blkLimits(LOW, JAXIS) .AND. & 
                                        k .eq. blkLimits(LOW, KAXIS)) then
                       temp_BoxVal(iter+1) = 0.0
                    end if

              !if (gr_meshMe .eq. 0 .AND. blockID .eq. 3 .AND. i .eq. 6 .AND. j.eq.10) then
              !   BoxVal(iter+1) = 999999.999999
              !end if

                 end if                 

              else  !- kpd - Boundary Node

                 MdensXL = 0.0
                 temp_BoxVal(iter+1) = 0.0

              end if
              
              !! i+1,j,k
              if (i /= blkLimits(HIGH, IAXIS) .or. (faces(2,IAXIS) == NOT_BOUNDARY)) then 
                 MdensXR = facevarx(MGW8_FACE_VAR,i+1,j  ,k)
                 !MdensXR = facevarx(RH1F_FACE_VAR,i+1,j  ,k) + facevarx(RH2F_FACE_VAR,i+1,j  ,k)
                 !print*,"Mdens is hardwired to 1 in createMatrix_KPD"
                 !MdensXR = 1.0 

                 if (MdensXR .le. 0.00 .OR. MdensXR .gt. 1300.0) then
                    print*,"ERROR in hypreCreateMatrix: Density is out of physical bounds!"
                 end if


                 if (i == blkLimits(HIGH, IAXIS) .and. mylevel /= gr_hypreNeghLevels(lb, RIGHT_EDGE, 1+K2D, 1+K3D)) then !! F/C boundary.                     
                 else                        
                    !! STENCILED RELATIONSHIP.                       

                    !*********************************************
                    !- kpd - For FLASH4 implementation for i+1,j,k
                    !*********************************************
                    temp_BoxVal(iter+2) = -1.d0 / (del(IAXIS)**2.0) * MdensXR
                    if (lb .eq. 1 .AND. gr_meshMe .eq. 0 .AND.             &
                                        i .eq. blkLimits(LOW, IAXIS) .AND. &
                                        j .eq. blkLimits(LOW, JAXIS) .AND. &
                                        k .eq. blkLimits(LOW, KAXIS)) then
                       temp_BoxVal(iter+2) = 0.0
                    end if
                 end if

              else  !- kpd - Boundary Node

                 MdensXR = 0.0
                 temp_BoxVal(iter+2) = 0.0

              end if
              
              !*******************************************
              !- kpd - For FLASH4 implementation for i,k,j
              !*******************************************
              temp_BoxVal(iter) = temp_BoxVal(iter) + 1./(del(IAXIS)**2.0) *(MdensXL + MdensXR)                             

              !print*,"AMAT",i,j,k,del(IAXIS),MdensXL,MdensXR,BoxVal(iter),BoxVal(iter+1),BoxVal(iter+2)
              
#if NDIM >= 2
              condjmh = 0.
              condjph = 0.
              
              !! i,j-1,k
              if ((j /= blkLimits(LOW, JAXIS)) .or. (faces(1,JAXIS) == NOT_BOUNDARY)) then                  

                 MdensYL = facevary(MGW8_FACE_VAR,i  ,j  ,k)
                 !MdensYL = facevary(RH1F_FACE_VAR,i  ,j  ,k)+facevary(RH2F_FACE_VAR,i  ,j  ,k)
                 !print*,"Mdens is hardwired to 1 in createMatrix_KPD"
                 !MdensYL = 1.0 

                 if (MdensYL .le. 0.00 .OR. MdensYL .gt. 1300.0) then
                    print*,"ERROR in hypreCreateMatrix: Density is out of physical bounds!"
                 end if

                 if (j == blkLimits(LOW, JAXIS) .and. mylevel /= gr_hypreNeghLevels(lb,1+K2D,LEFT_EDGE,1+K3D)) then  !! F/C boundary.                     
                    
                 else
                    !! STENCILED RELATIONSHIP.                       
                    !*********************************************
                    !- kpd - For FLASH4 implementation for i,j-1,k
                    !*********************************************
                    temp_BoxVal(iter+3) = -1.0 / (del(JAXIS)**2.0) * MdensYL
                    if (lb .eq. 1 .AND. gr_meshMe .eq. 0 .AND.             &
                                        i .eq. blkLimits(LOW, IAXIS) .AND. &
                                        j .eq. blkLimits(LOW, JAXIS) .AND. &
                                        k .eq. blkLimits(LOW, KAXIS)) then
                       temp_BoxVal(iter+3) = 0.0
                    end if
                 end if

              else   !- kpd - Boundary Node

                 MdensYL = 0.0
                 temp_BoxVal(iter+3) = 0.0

              end if
              
              !! i,j+1,k
              if ((j /= blkLimits(HIGH, JAXIS)) .or. (faces(2,JAXIS) == NOT_BOUNDARY)) then                  

                 MdensYR = facevary(MGW8_FACE_VAR,i  ,j+1,k)
                 !MdensYR = facevary(RH1F_FACE_VAR,i  ,j+1,k)+facevary(RH2F_FACE_VAR,i  ,j+1,k)
                 !print*,"Mdens is hardwired to 1 in createMatrix_KPD"
                 !MdensYR = 1.0 

                 if (MdensYR .le. 0.00 .OR. MdensYR .gt. 1300.0) then
                    print*,"ERROR in hypreCreateMatrix: Density is out of physical bounds!"
                 end if


                 if (j ==  blkLimits(HIGH, JAXIS) .and. mylevel /= gr_hypreNeghLevels(lb,1+K2D,RIGHT_EDGE,1+K3D)) then !! F/C boundary.  
                 else
                    !! STENCILED RELATIONSHIP.                       

                    !*********************************************
                    !- kpd - For FLASH4 implementation for i,j+1,j
                    !*********************************************
                    temp_BoxVal(iter+4) = -1.0 / (del(JAXIS)**2.0) * MdensYR
                    if (lb .eq. 1 .AND. gr_meshMe .eq. 0 .AND.             &
                                        i .eq. blkLimits(LOW, IAXIS) .AND. &
                                        j .eq. blkLimits(LOW, JAXIS) .AND. &
                                        k .eq. blkLimits(LOW, KAXIS)) then
                       temp_BoxVal(iter+4) = 0.0
                    end if
                 end if

              else  !- kpd - Boundary Node
 
                 MdensYR = 0.0
                 temp_BoxVal(iter+4) = 0.0

              end if
              
!print*,"Creating",i,j,facevarx(MGW8_FACE_VAR,i ,j,k),facevarx(MGW8_FACE_VAR,i+1,j,k), &
!        facevary(MGW8_FACE_VAR,i  ,j,k),facevary(MGW8_FACE_VAR,i  ,j+1,k)

              !*******************************************
              !- kpd - For FLASH4 implementation for i,j,k
              !*******************************************
              temp_BoxVal(iter) = temp_BoxVal(iter) + 1./(del(JAXIS)**2.0) * (MdensYL + MdensYR) 

              
#if NDIM == 3
              condkmh = 0.
              condkph = 0.              
              !! i,j,k-1
              if ((k /= blkLimits(LOW, KAXIS)) .or. (faces(1,KAXIS) == NOT_BOUNDARY)) then                  

                 MdensZL = facevarz(MGW8_FACE_VAR,i  ,j  ,k)
                 !MdensZL = facevarz(RH1F_FACE_VAR,i  ,j  ,k)+facevarz(RH2F_FACE_VAR,i  ,j  ,k)
                 !print*,"Mdens is hardwired to 1 in createMatrix_KPD"
                 !MdensZL = 1.0 

                 if (MdensZL .le. 0.00 .OR. MdensZL .gt. 1300.0) then
                    print*,"ERROR in hypreCreateMatrix: ZL Density is out of physical bounds!",MdensZL
                 end if


                 if (k == blkLimits(LOW, KAXIS) .and. mylevel /= gr_hypreNeghLevels(lb,1+K2D,1+K2D,LEFT_EDGE)) then  !! F/C boundary.                     
                    
                 else
                    !! STENCILED RELATIONSHIP.                       

                    !*********************************************
                    !- kpd - For FLASH4 implementation for i,j-1,k
                    !*********************************************
                    temp_BoxVal(iter+5) = -1.0 / (del(KAXIS)**2.0) * MdensZL
                    if (lb .eq. 1 .AND. gr_meshMe .eq. 0 .AND.             &
                                        i .eq. blkLimits(LOW, IAXIS) .AND. &
                                        j .eq. blkLimits(LOW, JAXIS) .AND. &
                                        k .eq. blkLimits(LOW, KAXIS)) then
                       temp_BoxVal(iter+5) = 0.0
                    end if
                 end if

              else  !- kpd - Boundary Node

                 MdensZL = 0.0
                 temp_BoxVal(iter+5) = 0.0

              end if
              
              !! i,j,k+1
              if ((k /= blkLimits(HIGH,KAXIS)) .or. (faces(2,KAXIS) == NOT_BOUNDARY)) then                  

                 MdensZR = facevarz(MGW8_FACE_VAR,i  ,j,k+1)
                 !MdensZR = facevarz(RH1F_FACE_VAR,i  ,j,k+1)+facevarz(RH2F_FACE_VAR,i  ,j,k+1)
                 !print*,"Mdens is hardwired to 1 in createMatrix_KPD"
                 !MdensZR = 1.0 

                 if (MdensZR .le. 0.00099 .OR. MdensZR .gt. 1300.0) then
                    print*,"ERROR in hypreCreateMatrix: ZR Density is out of physical bounds!",MdensZR
                 end if

                 if (k ==  blkLimits(HIGH,KAXIS) .and. mylevel /= gr_hypreNeghLevels(lb,1+K2D,1+K2D,RIGHT_EDGE)) then !! F/C boundary.                     
                 else
                    !! STENCILED RELATIONSHIP.                       

                    !*********************************************
                    !- kpd - For FLASH4 implementation for i,j-1,k
                    !*********************************************
                    temp_BoxVal(iter+6) = -1.0 / (del(KAXIS)**2.0) * MdensZR
                    if (lb .eq. 1 .AND. gr_meshMe .eq. 0 .AND.             &
                                        i .eq. blkLimits(LOW, IAXIS) .AND. &
                                        j .eq. blkLimits(LOW, JAXIS) .AND. &
                                        k .eq. blkLimits(LOW, KAXIS)) then
                       temp_BoxVal(iter+6) = 0.0
                    end if
                 end if

              else  !- kpd - Boundary Node

                 MdensZR = 0.0
                 temp_BoxVal(iter+6) = 0.0

              end if
              
              !*******************************************
              !- kpd - For FLASH4 implementation for i,j,k
              !*******************************************
              temp_BoxVal(iter) = temp_BoxVal(iter) + 1./(del(KAXIS)**2.0) * (MdensZL + MdensZR)

#endif              
#endif

              iter = iter + nentries

           end do
        end do
     end do     
     
!print*,"CREATE",gr_meshMe,blockID,gr_hypreLower(lb,1:NDIM),gr_hypreUpper(lb,1:NDIM),nentries,stencil_indices(1:nentries)
!do i=1,iter
!   print*,i,BoxVal(i),nentries,product(datasize(1:NDIM)),datasize(1),datasize(2),datasize(2)
!end do

!    if (NDIM .eq. 2) then
     !call HYPRE_SStructMatrixSetBoxValues(gr_hypreMatA, mypart, gr_hypreLower(lb,1:NDIM), &
     !     gr_hypreUpper(lb,1:NDIM), var, nentries, stencil_indices(1:nentries), BoxVal(:), ierr)
     call HYPRE_SStructMatrixSetBoxValues(gr_hypreMatA, mypart, gr_hypreLower(lb,1:NDIM), &
          gr_hypreUpper(lb,1:NDIM), var, nentries, stencil_indices(1:nentries), temp_BoxVal, ierr)
!    else
!     call HYPRE_SStructMatrixSetBoxValues(gr_hypreMatA, mypart, gr_hypreLower(lb,1:NDIM), &
!          gr_hypreUpper(lb,1:NDIM), var, nentries, stencil_indices(1:nentries), BoxVal(1:28672), ierr)
!    end if
      
     
    
     call Timers_start("gr_hypreApplyBcToEdge")     
     dir = ILO_FACE
     do i = IAXIS, NDIM      !- kpd - Go through x-dir, y-dir, & z-dir
        do j = LOW, HIGH     !        Go through 1,2
        !print*,"Boundary i,j",i,j,faces(j,i),NOT_BOUNDARY,LOW,HIGH
           !!- kpd - Is this a problem ??? Never goes in here ???
           !if (faces(j,i) /= NOT_BOUNDARY) then               
           !   call gr_hypreApplyBcToFace(blkLimits,blkLimitsGC,mypart,var,iFactorB,bcTypes(dir),dir, &
           !        bcValues(:,dir), dt, theta, del(i), gr_hypreLower(lb,:), dirichlet_multiplier, faceAreas(i,:,:,:), solnVec)
           !end if

           !   !- kpd - Gets called once per boundary (4 times in 2D)
           !   call gr_hypreApplyBcToFace_KPD(blkLimits,blkLimitsGC,mypart,var,iFactorB,bcTypes(dir),dir, &
           !        bcValues(:,dir), dt, theta, del(i), gr_hypreLower(lb,:), dirichlet_multiplier, faceAreas(i,:,:,:), solnVec)
           dir = dir + 1
        end do
     end do
     call Timers_stop("gr_hypreApplyBcToEdge") 
     
     call Grid_releaseBlkPtr(blockID, solnVec)             

     call Grid_releaseBlkPtr(blockID,facevarx ,FACEX)
     call Grid_releaseBlkPtr(blockID,facevary ,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facevarz ,FACEZ)
#endif

!     deallocate (faceAreas)
!     deallocate (BoxVal)
     
     end if  !- kpd - End level = lref if
  end do     !- kpd - End block loop
  !******************************************************************************************
  !******************************************************************************************
 

 !!-----------------------------------------------------------------------
 !!         THIS IS A GLOBAL CALL.
 !
 ! - kpd - Here the Amat is Asembled!
 !!-----------------------------------------------------------------------
 call HYPRE_SStructMatrixAssemble(gr_hypreMatA, ierr)
 
 if (blockCount > 0) then
!    deallocate (xflux)
!    deallocate (yflux)
!    deallocate (zflux)
 end if 
 
 call Timers_stop("gr_hypreCreateMatrix") 
 
 return
 
end subroutine gr_hypreCreateMatrix_KPD
