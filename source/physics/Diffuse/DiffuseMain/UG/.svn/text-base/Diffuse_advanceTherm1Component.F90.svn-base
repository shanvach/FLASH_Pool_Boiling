!!****if* source/physics/Diffuse/DiffuseMain/UG/Diffuse_advanceTherm1Component
!!
!!  NAME 
!!
!!  Diffuse_advanceTherm1Component
!!
!!  SYNOPSIS
!!
!!  call Diffuse_advanceTherm1Component(integer, intent(IN)           :: iVar,
!!                                      integer, intent(IN)           :: iFactorB,
!!                                      integer, intent(IN)           :: iFactorA,
!!                                      integer, intent(IN)           :: bcTypes(6),
!!                                      real,    intent(IN)           :: bcValues(2,6),
!!                                      real,    intent(IN)           :: dt,
!!                                      real,    intent(IN)           :: scaleFact,
!!                                      real,              intent(IN) :: chi,
!!                                      real,              intent(IN) :: theta,
!!                                      integer, OPTIONAL, intent(IN) :: pass,
!!                                      integer,           intent(IN) :: blockCount,
!!                                      integer,dimension(blockCount),intent(IN) :: blockList,
!!                                      integer, intent(IN), OPTIONAL :: iFactorC,
!!                                      integer, intent(IN), OPTIONAL :: iFactorD)
!!
!!
!!  DESCRIPTION 
!!      This routine advances the diffusion equation  of the form,
!!
!!      A * dV/dt = d/dx(B*dV/dx) + d/dy(B*dV/dy) + d/dx(B*dV/dz) + C*V + D
!!
!!
!!
!! ARGUMENTS
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!  
!!
!!***

!!REORDER(4): solnVec


subroutine Diffuse_advanceTherm1Component(iVar, iFactorB, iFactorA, bcTypes, bcValues, &
     dt, scaleFact, chi, theta, pass, blockCount, blockList, iFactorC, iFactorD)

 use diff_saData,       ONLY : diff_scaleFactThermSaTime
 use Diffuse_data,      ONLY : diff_meshMe,diff_geometry,diff_meshComm
 use Grid_interface,    ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
                               Grid_getBlkIndexLimits, Grid_fillGuardCells
 
 use Timers_interface,  ONLY : Timers_start, Timers_stop
 
 use diff_interface, ONLY    : diff_applyPC, diff_computeAAMatrix, &
                               diff_computeILU, diff_dotproduct, diff_computeB, &
                               diff_computeAX
 
 implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"
#include "Eos.h"

 integer,                      intent(IN) :: iVar
 integer,                      intent(IN) :: iFactorB
 integer,                      intent(IN) :: iFactorA
 integer, OPTIONAL,            intent(IN) :: iFactorC
 integer, OPTIONAL,            intent(IN) :: iFactorD
 integer,                      intent(IN) :: bcTypes(6)
 real,                         intent(IN) :: bcValues(2,6)
 real,                         intent(IN) :: dt
 real,                         intent(IN) :: scaleFact
 real,                         intent(IN) :: chi
 real,                         intent(IN) :: theta
 integer, OPTIONAL,            intent(IN) :: pass
 integer,                      intent(IN) :: blockCount
 integer,dimension(blockCount),intent(IN) :: blockList

 ! Local
 integer                   :: i, j, k, lb, ierr
 real, dimension(MDIM)     :: del 
 real                      :: CDiv 
 logical                   :: mask(NUNK_VARS)
 real, POINTER, DIMENSION(:,:,:,:) :: solnVec
 integer, dimension(2,MDIM)        :: blkLimitsGC, blkLimits 
 real    :: alpha, Beta, APdotP(blockCount), newRdotZ(blockCount), RdotZ(blockCount)
 integer :: datasize(MDIM), cgiter,datasizeGC(MDIM)
 real, allocatable, dimension(:,:,:,:):: AP
 real, allocatable, dimension(:,:,:,:):: R, Z 
 real    :: GlobalRdotZ, GlobalAPdotP,GlobalnewRdotZ


 !! CG Exit conditions, mimics ones in PetSc.
 real    :: diff_MaxIteration = 10000   ! Max allowed iterations.
 real    :: diff_atol         = 1.0E-16 ! Absolute tolerence.
 real    :: diff_rtol         = 1.0E-6  ! Relative tolerence.
 real    :: diff_dtol         = 1.0E+5  ! Divergence tolerence.
 logical :: converged
 real    :: InitialResidue
 
 logical :: diff_saveAAMat = .TRUE.
 logical :: diff_usePrecon = .TRUE.
 
 !! for CSR A Matrix.
 integer :: N, NZ
 integer, allocatable,dimension(:)   :: IA, JA, UPTR
 real, allocatable,dimension(:,:) :: AAs
 
 
 call Timers_start("Diffuse_advanceTherm1Component")
 
 call Grid_fillGuardCells(CENTER,ALLDIR) 
 
 ! GC needs to be exchanged for Krylov vector only, Mask the rest.
 mask           = .FALSE.
 mask(DTMP_VAR) = .TRUE. ! Krylov Vector 
 
 
 !!This is a minor hack, helps to get the max block size (typically 8X8X8)
 !!This might change in chombo and can be ignored for now.
 datasize  (1:MDIM) = 0.
 datasizeGC(1:MDIM) = 0.
 do lb = 1, blockCount           
    call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)                      
    datasizeGC(1:MDIM)=max(blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1,datasizeGC(1:MDIM))    
    datasize  (1:MDIM)=max(blkLimits  (HIGH,1:MDIM)-blkLimits  (LOW,1:MDIM)+1,datasize  (1:MDIM))    

 end do
 

 !! N  - Size of A matrix
 !! NZ - Number of non zero elements in A matrix. 
 N = PRODUCT(datasize)
#if NDIM == 1 
 NZ = 3*(datasize(IAXIS)-2) + 2*2
#else
 NZ = 5*(datasize(IAXIS)-2)*(datasize (JAXIS)-2) + 4*(2*(datasize(JAXIS)-2)+2*(datasize(IAXIS)-2)) + 3*4  
#endif
!! DEV: NDIM==3 case is missing.
 
 if (diff_usePrecon) then 
    allocate (JA(NZ))
    allocate (IA(N+1))
    allocate (UPTR(N))
    
    if (diff_saveAAMat) then
       allocate (AAs(NZ,blockCount))         
       !! Compute LU decomposition for each block and store it.    
       do lb = 1, blockCount    
          call diff_computeAAMatrix(AAs(:,lb), JA, IA, UPTR, NZ, N, blockList(lb), iFactorA, iFactorB, dt,theta,iFactorC) 
          call diff_computeILU(AAs(:,lb), JA, IA, UPTR, NZ, N)      
       end do
    else    
       !! Compute LU decomposition for each block as and when needed.
       allocate (AAs(NZ,1))     
    end if
 end if
    
 !! For Residue,
 allocate(R  (blockCount,datasizeGC(IAXIS),datasizeGC(JAXIS), datasizeGC(KAXIS)))
 
 !! For preconditioned systems, MZ=R, M ~ A, else Z = R
 allocate(Z  (blockCount,datasizeGC(IAXIS),datasizeGC(JAXIS), datasizeGC(KAXIS))) 
 
 !! Matrix vector product, A X P, here we don't store A matrix.
 allocate(AP (blockCount,datasizeGC(IAXIS),datasizeGC(JAXIS), datasizeGC(KAXIS))) 
 
 !! Compute RHS of AX = B
 call diff_computeB (blockCount, blockList, iVar, iFactorA, iFactorB, dt, theta, bcTypes, bcValues, iFactorD)

 
 R  = 0.
 AP = 0.
 Z  = 0.
 
 !! Start CG iterations.
 cgiter = 1
 
 !!===============================================================================================
 !! Obtain the first Krylov vector, 
 !! This steps gets CG rolling by computing initial Residue, Krylov Vector.
 !!=============================================================================================== 
 do lb = 1, blockCount         
    call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)    
    call Grid_getBlkPtr(blocklist(lb), solnVec)      

    call diff_computeAX (blockList(lb),iVar, AP(lb,:,:,:), iFactorB, iFactorA, blkLimits, blkLimitsGC, dt, theta,iFactorC)     
    
    
    R(lb,:,:,:) = SolnVec(WTMP_VAR,:,:,:)- AP(lb,:,:,:)                     
    
    Z(lb,:,:,:) = 0.0

    if (diff_usePrecon) then    
       if (diff_saveAAMat) then
          call diff_applyPC (Z(lb,:,:,:), R(lb,:,:,:), NZ, N, AAs(:,lb), JA, IA, UPTR, blkLimits, blkLimitsGC)
       else       
          call diff_computeAAMatrix(AAs(:,1), JA, IA, UPTR, NZ, N, blockList(lb), iFactorA, iFactorB, dt,theta,iFactorC)
          call diff_computeILU(AAs(:,1), JA, IA, UPTR, NZ, N)       
          call diff_applyPC (Z(lb,:,:,:), R(lb,:,:,:), NZ, N, AAs(:,1), JA, IA, UPTR, blkLimits, blkLimitsGC)
       end if
    else
       Z(lb,:,:,:) = R(lb,:,:,:)
    endif    

    solnVec(DTMP_VAR,:,:,:) = Z(lb,:,:,:)    ! P0 = Z0
    
    call diff_dotproduct(R(lb,:,:,:), Z(lb,:,:,:), RdotZ(lb), blkLimits, blkLimitsGC)     ! rj dot zj  
    
    call Grid_releaseBlkPtr(blocklist(lb), solnVec)
    
 end do !block loop
 !!===============================================================================================  

 converged = .FALSE.  
  
 !! Lets do a Global residue estimate, To start of the CG iterations.
 call mpi_allreduce (sum(RdotZ),  GlobalRdotZ, 1, FLASH_REAL, MPI_SUM, diff_meshComm, ierr) 
 call MPI_Barrier (diff_meshComm, ierr)
 
 InitialResidue = GlobalRdotZ
 
 if (sqrt(InitialResidue) > diff_atol) then          
    
    do while (.not. converged)       
       !! This GC fill helps deal with block/physical boundaries, only exchanges Krylov vector.    
       call Grid_fillGuardCells(CENTER,ALLDIR, masksize=NUNK_VARS, mask=mask)                 
      
       !!===============================================================================================
       !! Compute AP dot P across all blocks and processors.
       !!===============================================================================================
       APdotP        = 0.0
       GlobalAPdotP  = 0.0

       do lb = 1, blockCount     
          call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)    
          call Grid_getBlkPtr(blocklist(lb), solnVec)        
          
          !! compute Matrix vector product, AP         
          call diff_computeAX (blockList(lb),DTMP_VAR, AP(lb,:,:,:), iFactorB, iFactorA, blkLimits, blkLimitsGC, dt, theta,iFactorC)                    
          
          !! compute inner product (P,AP), stored in APdotP
          call diff_dotproduct(AP(lb,:,:,:), solnVec(DTMP_VAR,:,:,:), APdotP(lb),blkLimits, blkLimitsGC) 
          
          call Grid_releaseBlkPtr(blocklist(lb), solnVec)        
       end do
       
       !! Global reduction of APdotP.
       call mpi_allreduce (sum(APdotP), GlobalAPdotP, 1, FLASH_REAL, MPI_SUM, diff_meshComm, ierr)    
       call MPI_Barrier (diff_meshComm, ierr)
       !!===============================================================================================      
       
       alpha = GlobalRdotZ/GlobalAPdotP            
       
       !!===============================================================================================       
       !! Now with the new global alpha we 
       !! 1) update Solution (in unk). 
       !! 2) Compute updated residue per block.
       !! 3) Compute new R dot Z.
       !!===============================================================================================       
       do lb = 1, blockCount       
          call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)           
          call Grid_getBlkPtr(blocklist(lb), solnVec)           
          
          solnVec(iVar,:,:,:) =  solnVec(iVar,:,:,:)  + alpha*solnVec(DTMP_VAR,:,:,:)                  
          R(lb,:,:,:) = R(lb,:,:,:) - alpha*AP(lb,:,:,:)
          
          !! Now we need to compute the new 'P' vector.                    
          if (diff_usePrecon) then
             if (diff_saveAAMat) then
                call diff_applyPC (Z(lb,:,:,:), R(lb,:,:,:), NZ, N, AAs(:,lb), JA, IA, UPTR, blkLimits, blkLimitsGC)
             else       
                call diff_computeAAMatrix(AAs(:,1), JA, IA, UPTR, NZ, N, blockList(lb), iFactorA, iFactorB, dt, theta,iFactorC)
                call diff_computeILU(AAs(:,1), JA, IA, UPTR, NZ,N)       
                call diff_applyPC (Z(lb,:,:,:), R(lb,:,:,:), NZ, N, AAs(:,1), JA, IA, UPTR, blkLimits, blkLimitsGC)
             end if
          else
             Z(lb,:,:,:) = R(lb,:,:,:)
          endif
          
          !!Inner product of newly computed residue.
          call diff_dotproduct(R(lb,:,:,:),Z(lb,:,:,:),newRdotZ(lb),blkLimits, blkLimitsGC)       
          
          call Grid_releaseBlkPtr(blocklist(lb), solnVec)
       end do
       
       call mpi_allreduce (sum(newRdotZ), GlobalnewRdotZ, 1, FLASH_REAL, MPI_SUM, diff_meshComm, ierr)   
       call MPI_Barrier (diff_meshComm, ierr)
       !!===============================================================================================    
       
       Beta = GlobalnewRdotZ / GlobalRdotZ                 
       
       !!===============================================================================================       
       !!Update Krylov vector, i.e the next search direction.
       !!===============================================================================================       
       do lb = 1, blockCount       
          call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)           
          call Grid_getBlkPtr(blocklist(lb), solnVec)                                     
          solnVec(DTMP_VAR,:,:,:) = Z(lb,:,:,:) + Beta*solnVec(DTMP_VAR,:,:,:)          
          call Grid_releaseBlkPtr(blocklist(lb), solnVec)    
       end do
       !!===============================================================================================              
       
       !!Update old residue with new residue, increase number of iterations.       
       GlobalRdotZ = GlobalnewRdotZ        
       cgiter = cgiter + 1       

       !! check for convergence, diverenge or Max iterations.
       if (sqrt(GlobalRdotZ) <= MAX(diff_rtol*sqrt(initialResidue), diff_atol) .or. cgiter > diff_MaxIteration) then
          converged = .TRUE.  
       else
          if (sqrt(GlobalRdotZ) >= diff_dtol*sqrt(initialResidue)) then                          
             call Driver_abortFlash("Conjugate Gradient residue is diverging!") 
          endif
       end if
    end do ! CG iterations.    
    

    !if ( diff_meshMe == MASTER_PE ) then
    !   write(*,*) "CG Iterations:", cgiter, sqrt(GlobalRdotZ/InitialResidue), sqrt(InitialResidue)
    !end if
    
 end if
 
 
 !!Check for Max iterations, abort if more then set number of iterations.
 if ( diff_meshMe == MASTER_PE ) then
    if (cgiter > diff_MaxIteration) then
       call Driver_abortFlash("Conjugate Gradient: No convergence after Max Iterations, try higher value!") 
    end if
 end if 
 
 !! Flooring, avoids small (<< 0.0) negative floating point numbers.
 do lb = 1, blockCount
    call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blocklist(lb), solnVec)    
    solnVec(iVar,:,:,:) = max(solnVec(iVar,:,:,:), 1.0E-10) 
    call Grid_releaseBlkPtr(blocklist(lb), solnVec)
 end do
 
 !!Deallocate memory across the block.
 deallocate (R)
 deallocate (AP)
 deallocate (Z)
 
 if (diff_usePrecon) then
    deallocate (JA)
    deallocate (IA)
    deallocate (UPTR)
    deallocate (AAs) 
 end if
 
 call Timers_stop("Diffuse_advanceTherm1Component") 
 
 return
 
end subroutine Diffuse_advanceTherm1Component
