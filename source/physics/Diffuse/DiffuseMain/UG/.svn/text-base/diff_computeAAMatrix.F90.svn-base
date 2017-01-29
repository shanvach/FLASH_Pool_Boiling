!!****if* source/physics/Diffuse/DiffuseMain/UG/diff_computeAAMatrix
!!
!!  NAME 
!!
!!  diff_computeAAMatrix
!!
!!  SYNOPSIS
!!
!!  call diff_computeAAMatrix(AA, JA, IA, UPTR, NZ, N, blockID, iFactorA, iFactorB, iFactorC)!!
!!
!!  DESCRIPTION 
!!
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

subroutine diff_computeAAMatrix(AA, JA, IA, UPTR, NZ, N, blockID, iFactorA, iFactorB, dt,theta, iFactorC)  
  
  use Grid_interface,    ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
                                Grid_getBlkIndexLimits, &
                                Grid_getDeltas, Grid_getBlkBC
  use Timers_interface,  ONLY : Timers_start, Timers_stop
  
  implicit none
  
#include "Flash.h"
#include "constants.h"  
  
  integer, intent (IN)  :: blockID
  integer, intent (IN)  :: NZ
  integer, intent (IN)  :: N 
  real,    intent (IN)  :: dt
  real,    intent (IN)  :: theta  
  real,    intent (OUT) :: AA  (NZ)
  integer,    intent (OUT) :: JA  (NZ)
  integer,    intent (OUT) :: IA  (N+1)
  integer,    intent (OUT) :: UPTR(N)
  
  integer, intent(IN)           :: iFactorA
  integer, intent(IN)           :: iFactorB
  integer, OPTIONAL, intent(IN) :: iFactorC  
  
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer, DIMENSION(2,MDIM) :: blkLimitsGC, blkLimits 
  real, dimension(MDIM) :: del   

  real    :: CDiv
  integer :: pos_ijk, ia_iter, cs_iter
    
  real   :: condiph, condimh
  real   :: condjph, condjmh
  real   :: condkph, condkmh 

  integer :: i, j, k, pos
  integer :: datasize(MDIM), faces(2,MDIM)
  
  real    :: add2diag 
  
  call Timers_start("diff_computeAAMatrix")  
  
  cs_iter = 1
  ia_iter = 1
  
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)    
  call Grid_getDeltas(blockID, del)
  call Grid_getBlkPtr(blockID, solnVec)      
  call Grid_getBlkBC(blockID, faces)

  datasize  (1:MDIM)=blkLimits  (HIGH,1:MDIM)-blkLimits  (LOW,1:MDIM)+1   
  
  do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
     do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)            
           
           cDiv   = 1.0/solnVec(iFactorA,i,j,k)   
           
           add2diag = 0.0
           
           !! i,j,k => position in matrix            
           pos_ijk = (i-blkLimits(LOW,IAXIS)+1) + datasize(IAXIS)*(j-1-blkLimits(LOW,JAXIS)+1) + &
                datasize(IAXIS)*datasize(JAXIS)*(k-1-blkLimits(LOW,KAXIS)+1)
           
           condiph = 0.5*(solnVec(iFactorB,i+1,j,k)+ solnVec(iFactorB,i,j,k))*cDiv
           condimh = 0.5*(solnVec(iFactorB,i-1,j,k)+ solnVec(iFactorB,i,j,k))*cDiv                                              
#if NDIM >= 2             
           condjph = 0.5*(solnVec(iFactorB,i,j+1,k)+ solnVec(iFactorB,i,j,k))*cDiv
           condjmh = 0.5*(solnVec(iFactorB,i,j-1,k)+ solnVec(iFactorB,i,j,k))*cDiv             
#if NDIM == 3             
           condkph = 0.5*(solnVec(iFactorB,i,j,k+1)+ solnVec(iFactorB,i,j,k))*cDiv
           condkmh = 0.5*(solnVec(iFactorB,i,j,k-1)+ solnVec(iFactorB,i,j,k))*cDiv 
#endif  
#endif
           IA(ia_iter) = cs_iter    
           
#if NDIM == 3
           !! i,j,k-1
           if (k /= blkLimits(LOW, KAXIS)) then              
              AA(cs_iter) = -condkmh*theta*dt/(del(3)**2)
              pos = pos_ijk - datasize(IAXIS)*datasize(JAXIS)
              JA(cs_iter) = pos
              cs_iter = cs_iter + 1                            
           else
              if (faces(1,KAXIS) == OUTFLOW .or. faces(1,KAXIS) == REFLECTING) then                 
                 add2diag = add2diag -condkmh*theta*dt/(del(3)**2)              
              end if
           end if
#endif
           
#if NDIM == 2           
           !! i,j-1,k
           if (j /= blkLimits(LOW, JAXIS)) then              
              AA(cs_iter) = -condjmh*theta*dt/(del(2)**2)
              pos = pos_ijk - datasize(IAXIS)
              JA(cs_iter) = pos
              cs_iter = cs_iter + 1
           else
              if (faces(1,JAXIS) == OUTFLOW .or. faces(1,JAXIS) == REFLECTING) then                 
                 add2diag = add2diag -condjmh*theta*dt/(del(2)**2)              
              end if
           end if
#endif
           !! i-1,j,k
           if (i /= blkLimits(LOW, IAXIS)) then              
              AA(cs_iter) = -condimh*theta*dt/(del(1)**2)
              pos = pos_ijk - 1
              JA(cs_iter) = pos
              cs_iter = cs_iter + 1               
           else
              if (faces(1,IAXIS) == OUTFLOW .or. faces(1,IAXIS) == REFLECTING) then
                 add2diag = add2diag - condimh*theta*dt/(del(1)**2)
              end if
           end if
           
           !! i,j,k
           AA(cs_iter) = 1.0 + theta*(Condimh+Condiph)*dt/(del(1)*del(1))
#if NDIM >= 2                
           AA(cs_iter) =  AA(cs_iter) + theta*(Condjmh+Condjph)*dt/(del(2)*del(2))
#if NDIM == 3                         
           AA(cs_iter) =  AA(cs_iter) + theta*(Condkmh+Condkph)*dt/(del(3)*del(3))               
#endif  
#endif  
           if (present(iFactorC)) then               
              AA(cs_iter) = AA(cs_iter) + dt * solnVec(iFactorC,i,j,k)
           end if
           
           UPTR(ia_iter) = cs_iter
           pos = pos_ijk
           JA(cs_iter) = pos            
           cs_iter = cs_iter + 1           
           
           !! i+1,j,k
           if (i /= blkLimits(HIGH, IAXIS)) then              
              AA(cs_iter) = -condiph*theta*dt/(del(1)*del(1))  
              pos = pos_ijk + 1
              JA(cs_iter) = pos
              cs_iter = cs_iter + 1   
           else
              if (faces(2,IAXIS) == OUTFLOW .or. faces(2,IAXIS) == REFLECTING) then
                 add2diag = add2diag -condiph*theta*dt/(del(1)*del(1))
              end if
           end if
           
#if NDIM == 2
           !! i,j+1,k
           if (j /= blkLimits(HIGH, JAXIS)) then              
              AA(cs_iter) = -condjph*theta*dt/(del(2)*del(2))  
              pos = pos_ijk + datasize(IAXIS)
              JA(cs_iter) = pos
              cs_iter = cs_iter + 1  
           else
              if (faces(2,JAXIS) == OUTFLOW .or. faces(2,JAXIS) == REFLECTING) then
                 add2diag = add2diag -condjph*theta*dt/(del(2)*del(2)) 
              end if
           end if
#endif
           
           
#if NDIM == 3
           !! i,j,k+1
           if (k /= blkLimits(HIGH, KAXIS)) then
              AA(cs_iter) = -condkph*theta*dt/(del(3)*del(3))  
              pos = pos_ijk + datasize(IAXIS)*datasize(JAXIS)
              JA(cs_iter) = pos
              cs_iter = cs_iter + 1    
           else
              if (faces(2,KAXIS) == OUTFLOW .or. faces(2,KAXIS) == REFLECTING) then
                 add2diag = add2diag -condkph*theta*dt/(del(3)**2)              
              end if
           end if
#endif


           AA(UPTR(ia_iter)) =  AA(UPTR(ia_iter)) + add2diag
           
           ia_iter = ia_iter + 1
           
        end do
     end do
  end do
  
  IA(ia_iter) = 1 + NZ
  
  call Grid_releaseBlkPtr(blockID, solnVec)    

  
  call Timers_stop("diff_computeAAMatrix")    
  
end subroutine diff_computeAAMatrix
