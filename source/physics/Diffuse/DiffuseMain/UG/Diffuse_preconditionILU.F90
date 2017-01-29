!!****if* source/physics/Diffuse/DiffuseMain/UG/Diffuse_preconditionILU
!!
!!  NAME 
!!
!!  Diffuse_advanceTherm1Component
!!
!!  SYNOPSIS
!!
!!  call Diffuse_advanceTherm1Component(integer, intent(IN)           :: iVar                                              
!!                                      integer, intent(IN)           :: iFactorB
!!                                      integer, intent(IN)           :: iFactorA
!!                                      integer, intent(IN)           :: bcTypes(6)
!!                                      real,    intent(IN)           :: bcValues(2,6)
!!                                      real,    intent(IN)           :: dt
!!                                      real,    intent(IN)           :: scaleFact
!!                                      real,    OPTIONAL, intent(IN) :: chi
!!                                      integer, OPTIONAL, intent(IN) :: pass
!!                                      integer,           intent(IN) :: blockCount
!!                                      integer,dimension(blockCount),intent(IN) :: blockList
!!                                      integer, intent(IN), OPTIONAL :: iFactorC
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
!! NOTES:
!!  
!!
!!***

!!REORDER(4): solnVec

subroutine Diffuse_precondition(Z, Res, blkLimits, blkLimitsGC, iFactorB, iFactorA, del, dt, theta,iFactorC)

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
  
  real, intent(OUT)    :: Z  (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                              blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))  
  real, intent(IN)     :: Res(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))  
  real, intent(IN)     :: iFactorB (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) 
  real, intent(IN)     :: iFactorA (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
  real, intent(IN), OPTIONAL     :: iFactorC (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))

  real, dimension(MDIM), intent(IN)            :: del
  real, intent(IN)                             :: dt
  real, intent (IN)                            :: theta  
  
  real, allocatable,dimension(:) :: LD,DD,UD,S,X
  integer :: datasize(MDIM)
  integer :: i, j, k
  real    :: Cond_R, Cond_L,CDiv
  
  real, allocatable,dimension(:,:) :: M
  
  
  real   :: condiph, condimh
  real   :: condjph, condjmh
  real   :: condkph, condkmh
  
  integer :: N, NZ, cs_iter, ia_iter
  real, allocatable,dimension(:) :: IA, JA, UPTR, AA, LU, work
  integer :: pos, pos_ijk,k1,k2,jj
  logical :: diagreached
  

  Z = 0.0
  
  
#ifdef ILU
  
  !! The idea is to store the A matrix in CSR, then perform ILU on that matrix.
  !! This would be ILU(0) since it operates only on specific positions on CSR.
  
  datasize  (1:MDIM)=blkLimits  (HIGH,1:MDIM)-blkLimits  (LOW,1:MDIM)+1  
  
  !! Total number of interior cells.
  N = PRODUCT(datasize)
  
  !-----------------------------------------------
  !1D
  ! (datasize(IAXIS)-2) => 3 point stencil
  ! 2                   => 2 point stencil (corners)
  !2D
  !  (datasize(IAXIS)-2)*(datasize (JAXIS)-2) => 5 point stencil
  !2*(datasize(JAXIS)-2)                      => 4 point stencil 
  !2*(datasize(IAXIS)-2)                      => 4 point stencil
  !4                                          => 3 point stencil  
  !
  
  !! TODO: 3D
#if NDIM == 1 
  NZ = 3*(datasize(IAXIS)-2) + 2*2
#else
  NZ = 5*(datasize(IAXIS)-2)*(datasize (JAXIS)-2) + 4*(2*(datasize(JAXIS)-2)+2*(datasize(IAXIS)-2)) + 3*4  
#endif 

  allocate (AA(NZ))
  allocate (LU(NZ))
  allocate (JA(NZ))
  allocate (IA(N+1))
  allocate (UPTR(N))
  allocate (work(N))
  allocate (S(N))
  allocate (X(N))
  
  
!!$  !! Test problem 1
!!$  AA   = (/2.0,1,1,2.0,1,1,2.0,1,1,2.0,1,1,2.0/)  !! 13
!!$  JA   = (/1,2, 1,2,3,2,3,4,3,4,5, 4,5/)  !! 13
!!$  IA   = (/1,3,6,9,12,14/)              !! 6
!!$  uptr = (/1,4,7,10,13/)                !! 5
!!$
!!$
!!$  !! Test problem 2
!!$  AA   = (/5,1,1,1,1, &
!!$       1,5,1,1,1, &
!!$       1,1,5,1,1, &
!!$       1,1,1,5,1, &
!!$       1,1,1,1,5/)  
!!$  JA   = (/1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5/)
!!$  IA   = (/1,6,11,16,21,26/)
!!$  uptr = (/1,7,13,19,25/)          
  
  cs_iter = 1
  ia_iter = 1
  
  do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
     do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)            
           
           cDiv   = 1.0/iFactorA(i,j,k)                 
           
           !! i,j,k => position in matrix            
           pos_ijk = (i-blkLimits(LOW,IAXIS)+1) + datasize(IAXIS)*(j-1-blkLimits(LOW,JAXIS)+1) + &
                datasize(IAXIS)*datasize(JAXIS)*(k-1-blkLimits(LOW,KAXIS)+1)
           
           condiph = 0.5*(iFactorB(i+1,j,k)+ iFactorB(i,j,k))*cDiv
           condimh = 0.5*(iFactorB(i-1,j,k)+ iFactorB(i,j,k))*cDiv                                              
#if NDIM >= 2             
           condjph = 0.5*(iFactorB(i,j+1,k)+ iFactorB(i,j,k))*cDiv
           condjmh = 0.5*(iFactorB(i,j-1,k)+ iFactorB(i,j,k))*cDiv             
#if NDIM == 3             
           condkph = 0.5*(iFactorB(i,j,k+1)+ iFactorB(i,j,k))*cDiv
           condkmh = 0.5*(iFactorB(i,j,k-1)+ iFactorB(i,j,k))*cDiv 
#endif  
#endif
           IA(ia_iter) = cs_iter    
           
           !! i,j,k-1
           if (k /= blkLimits(LOW, KAXIS)) then              
              AA(cs_iter) = -condkmh*theta*dt/(del(3)**2)
              pos = pos_ijk - datasize(IAXIS)*datasize(JAXIS)
              JA(cs_iter) = pos
              cs_iter = cs_iter + 1                            
           end if
           
           !! i,j-1,k
           if (j /= blkLimits(LOW, JAXIS)) then              
              AA(cs_iter) = -condjmh*theta*dt/(del(2)**2)
              pos = pos_ijk - datasize(IAXIS)
              JA(cs_iter) = pos
              cs_iter = cs_iter + 1
           end if
            
           !! i-1,j,k
           if (i /= blkLimits(LOW, IAXIS)) then              
              AA(cs_iter) = -condimh*theta*dt/(del(1)**2)
              pos = pos_ijk - 1
              JA(cs_iter) = pos
              cs_iter = cs_iter + 1               
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
              AA(cs_iter) = AA(cs_iter) + dt * iFactorC (i,j,k)
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
           end if
    
           
           !! i,j+1,k
           if (j /= blkLimits(HIGH, JAXIS)) then              
              AA(cs_iter) = -condjph*theta*dt/(del(2)*del(2))  
              pos = pos_ijk + datasize(IAXIS)
              JA(cs_iter) = pos
              cs_iter = cs_iter + 1            
           end if
           
           !! i,j,k+1
           if (k /= blkLimits(HIGH, KAXIS)) then
              AA(cs_iter) = -condkph*theta*dt/(del(3)*del(3))  
              pos = pos_ijk + datasize(IAXIS)*datasize(JAXIS)
              JA(cs_iter) = pos
              cs_iter = cs_iter + 1            
           end if
           
           S (pos_ijk) = Res(i,j,k)                   
           
           ia_iter = ia_iter + 1
           
        end do
     end do
  end do
  
  IA(ia_iter) = 1 + NZ
  
  !! DO ILU factorization.  
  !! LU maintains the same sparsity of AA hence ILU(0)
  LU = AA 
  work = 0
  do i = 2, N           
     k1 = IA(i)
     k2 = IA(i+1)-1     
     diagreached = .false.

     do j=k1, k2
        work(JA(j)) = j        
     end do     

     do while (JA(k1) .lt. i .and. .not.(diagreached))                 
        LU(K1) = LU(K1) / LU(uptr(JA(K1)))                  
        
        if (LU(uptr(JA(K1))) == 0.0) then
           write(*,*) "zero pivot"
           pause
        end if
        
        do jj = uptr(JA(K1))+1, ia(JA(K1)+1)-1                       
           if (work(ja(jj)) /= 0) then 
              LU(work(JA(jj))) =  LU(work(JA(jj))) - LU(K1)*LU(jj)                                   
           end  if
        end do
        
        if (JA(k1) .eq.uptr(i)-1) then            
           diagreached = .true.           
        else                   
           k1 = k1 + 1           
        end if
        
     end do
     work = 0     
  end do


  
  !! with LU decomposition, we need to back/forward Solve.
  !! AX = B
  !! LUX = B, UX = Y
  !! LY = B -> Lower diagnol system.
  !! UX = Y -> Upper diagnol system.   
  X = 0.0
  
  ! LY = B  
  X(1) = S(1)
  do i = 2, N
     k1 = IA(i)
     k2 = IA(i+1)-1
     X(i)=S(i)- dot_product(LU(k1:k2),X(JA(k1:k2)))
  end do
  
  
  S(:) = X(:)
  ! UX = Y
  
  X(:) = 0.0

  X(N) = S(N) / LU(NZ)   
  
  do i=N-1, 1, -1
     k1 = IA(i+1)-1
     k2 = IA(i)
     X(i)=(S(i)- dot_product(LU(k2:k1),X(JA(k2:k1))))/LU(UPTR(i))
  end  do 


  !! Update solution.
  do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
     do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
           !! i,j,k => position in matrix            
           pos_ijk = (i-blkLimits(LOW,IAXIS)+1) + datasize(IAXIS)*(j-1-blkLimits(LOW,JAXIS)+1) + &
                datasize(IAXIS)*datasize(JAXIS)*(k-1-blkLimits(LOW,KAXIS)+1)
           
           Z(i,j,k) = X(pos_ijk)
           
        end do
     end do
  end do
  
  
  
  
  deallocate (AA)
  deallocate (LU)
  deallocate (JA)
  deallocate (IA)
  deallocate (UPTR)
  deallocate (S)
  deallocate (X)      
  deallocate (work)
  
#endif 
  
end subroutine Diffuse_precondition
