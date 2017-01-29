!!****if* source/physics/Diffuse/DiffuseMain/UG/Diffuse_precondition
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
  
  !! We need to solve MZ=R, where M ~ A   
  !#define USEGC

  datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1       
  
  allocate (LD(datasize(IAXIS)))  ! Lower diagonal.
  allocate (DD(datasize(IAXIS)))  ! Upper diagonal.
  allocate (UD(datasize(IAXIS)))  ! digaonal.     
  allocate ( S(datasize(IAXIS)))  ! RHS: Source Term.   
  allocate ( X(datasize(IAXIS)))    

!#define PC_INCOMPLETE_CHOLESKY
#define PC_THOMAS_1D
!#define PC_DIRECT_SOLVE_UG
!#define PC_BLOCK_JACOBI

!!-------------------------------------------------------------------------------------------------
#ifdef PC_THOMAS_1D
  do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
     cDiv   = 1.0/iFactorA(i,1,1)     
     Cond_L = 0.5*(iFactorB(i-1,1,1)+ iFactorB(i,1,1))*cDiv
     Cond_R = 0.5*(iFactorB(i+1,1,1)+ iFactorB(i,1,1))*cDiv       

     LD(i) = -Cond_L*theta*dt/(del(1)*del(1))     
     UD(i) = -Cond_R*theta*dt/(del(1)*del(1))         
     
     DD(i) = 1.0 + theta*(Cond_L+Cond_R)*dt/(del(1)*del(1))
     
     if (present(iFactorC)) then
        DD(i) = DD(i) + dt * iFactorC (i,1,1)
     end if
     
     S (i) = Res(i,1,1)       
  end do
  
  
!!$  !! dirichlet bc.
!!$  i = blkLimits(LOW, IAXIS)
!!$  S (i) = Res(i,1,1) - LD(i)*Z(i-1,1,1)
!!$
!!$  i = blkLimits(HIGH, IAXIS)
!!$  S (i) = Res(i,1,1) - UD(i)*Z(i+1,1,1)
  

  call Thomas3D(LD,DD,UD,S,blkLimits,blkLimitsGC)       
  
  do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
     Z(i,1,1) = S(i)
  end do

  
#endif
!!-------------------------------------------------------------------------------------------------  
#ifdef PC_DIRECT_SOLVE_UG
  do i = blkLimits(LOW, IAXIS)+1, blkLimits(HIGH, IAXIS)-1
     cDiv   = 1.0/iFactorA(i,1,1)
     
     Cond_L = 0.5*(iFactorB(i-1,1,1)+ iFactorB(i,1,1))*cDiv
     Cond_R = 0.5*(iFactorB(i+1,1,1)+ iFactorB(i,1,1))*cDiv       
     
     LD(i) = -Cond_L*theta*dt/(del(1)*del(1))
     
     UD(i) = -Cond_R*theta*dt/(del(1)*del(1))          
     
     DD(i) = (1.0 + theta*(Cond_L+Cond_R)*dt/(del(1)*del(1)))   
     
     if (present(iFactorC)) then
        DD(i) = DD(i) + dt * iFactorC (i,1,1)
     end if
     
     S (i) = Res(i,1,1)       
  end do
  
  ! Outflow
  i = blkLimits(LOW, IAXIS)
  cDiv   = 1.0/iFactorA(i,1,1)
  Cond_L = 0.5*(iFactorB(i-1,1,1)+ iFactorB(i,1,1))*cDiv
  Cond_R = 0.5*(iFactorB(i+1,1,1)+ iFactorB(i,1,1))*cDiv   
  
  
  LD(i) = -Cond_L*theta*dt/(del(1)*del(1))
  UD(i) = -Cond_R*theta*dt/(del(1)*del(1))                  
  DD(i) = (1.0 + theta*(Cond_L+Cond_R)*dt/(del(1)*del(1))) + LD(i)
  if (present(iFactorC)) then
     DD(i) = DD(i) + dt * iFactorC (i,1,1)
  end if
  S (i) = Res(i,1,1)

  ! Outflow
  i = blkLimits(HIGH, IAXIS)
  cDiv   = 1.0/iFactorA(i,1,1)
  Cond_L = 0.5*(iFactorB(i-1,1,1)+ iFactorB(i,1,1))*cDiv
  Cond_R = 0.5*(iFactorB(i+1,1,1)+ iFactorB(i,1,1))*cDiv
     
  LD(i) = -Cond_L*theta*dt/(del(1)*del(1))
  UD(i) = -Cond_R*theta*dt/(del(1)*del(1))     
  
  DD(i) = (1.0 + theta*(Cond_L+Cond_R)*dt/(del(1)*del(1))) + UD(i)
  if (present(iFactorC)) then
     DD(i) = DD(i) + dt * iFactorC (i,1,1)
  end if
  
  S (i) = Res(i,1,1)      
  
  call Thomas3D(LD,DD,UD,S,blkLimits,blkLimitsGC)     
  
  
  do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
     Z(i,1,1) = S(i)
  end do

#endif
!!-------------------------------------------------------------------------------------------------
#ifdef PC_INCOMPLETE_CHOLESKY
  allocate (M(datasize(IAXIS), datasize(IAXIS)))    

  M(:,:) = 0.0  

  do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
     cDiv   = 1.0/iFactorA(i,1,1)     
     Cond_L = 0.5*(iFactorB(i-1,1,1)+ iFactorB(i,1,1))*cDiv
     Cond_R = 0.5*(iFactorB(i+1,1,1)+ iFactorB(i,1,1))*cDiv        
     
     M(i,i  ) = 1.0 + theta*(Cond_L+Cond_R)*dt/(del(1)*del(1))
     M(i,i-1) = -Cond_L*theta*dt/(del(1)*del(1))
     M(i,i+1) = -Cond_R*theta*dt/(del(1)*del(1)) 

     if (present(iFactorC)) then
        M(i,i) = M(i,i) + dt * iFactorC (i,1,1)
     end if    
     
     S (i) = Res(i,1,1)           
     
  end do

  i = blkLimits(LOW, IAXIS)  
  S(i) = Res(i,1,1) - M(i,i-1)*Z(i-1,1,1)  
  

  i = blkLimits(HIGH, IAXIS)
  S(i) = Res(i,1,1) - M(i,i+1)*Z(i+1,1,1)
  
  call IncompleteCholesky(M(blkLimits(LOW, IAXIS):blkLimits(HIGH, IAXIS), &
                            blkLimits(LOW, IAXIS):blkLimits(HIGH, IAXIS)), &
                            blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1)

  

  !! change S(i) to deal with block boundaries.

  
  !! LUX = B, LY = B, UX = Y  
  !! Now we back/forward solve to obtain Z.

  z = 0.0
  !! 1) forward substitution. LY = B
  Z(blkLimits(LOW, IAXIS),1,1) =  S(blkLimits(LOW, IAXIS)) / M(blkLimits(LOW, IAXIS),blkLimits(LOW, IAXIS))
  do i = blkLimits(LOW, IAXIS)+1,  blkLimits(HIGH,IAXIS)   
     Z(i,1,1) = (S(i) - dot_product(M(i,i-1:i), Z(i-1:i,1,1)))/M(i,i)
  end do  
  S(:) = Z(:,1,1) 
  
  ! 2) back substitution, UX = Y
  Z(blkLimits(HIGH, IAXIS),1,1) = S(blkLimits(HIGH, IAXIS)) / M(blkLimits(HIGH, IAXIS),blkLimits(HIGH, IAXIS))
  do i = blkLimits(HIGH, IAXIS)-1,  blkLimits(LOW,IAXIS), -1
     Z(i,1,1) = (S(i) - dot_product(M(i+1:blkLimits(HIGH, IAXIS),i), Z(i+1:blkLimits(HIGH, IAXIS),1,1)))/M(i,i)
  end do  
  
  deallocate(M)
#endif
!!-------------------------------------------------------------------------------------------------
#ifdef PC_BLOCK_JACOBI

  do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
     cDiv   = 1.0/iFactorA(i,1,1)     
     Cond_L = 0.5*(iFactorB(i-1,1,1)+ iFactorB(i,1,1))*cDiv
     Cond_R = 0.5*(iFactorB(i+1,1,1)+ iFactorB(i,1,1))*cDiv        
     
     Diag = 1.0 + theta*(Cond_L+Cond_R)*dt/(del(1)*del(1))

     Z(i,1,1) = Res(i,1,1) / Diag          
  end do

#endif
!!-------------------------------------------------------------------------------------------------
  
  
  deallocate (LD)
  deallocate (UD)
  deallocate (DD)
  deallocate (S)
  deallocate (X)
  
end subroutine Diffuse_precondition


subroutine Thomas3D(UD,DD,LD,SS,blkLimits,blkLimitsGC)

  implicit none
  
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC

  real   , intent (IN)    :: LD (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS))
  real   , intent (IN)    :: UD (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS))
  real   , intent (INOUT) :: DD (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS))
  real   , intent (INOUT) :: SS (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS))

  real    :: R
  integer :: i,j


#ifdef USEGC
  do  i  = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)+1
     R     = UD(i)/DD(i-1)
     DD(i) = DD(i)-R*LD(i-1)
     SS(i) = SS(i)-R*SS(i-1)
  end do
  
  ! Back Substitution.
  SS (blkLimits(HIGH,IAXIS)+1) = SS(blkLimits(HIGH,IAXIS)+1)/DD(blkLimits(HIGH,IAXIS)+1)   
 
  do i    =  1, (blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS))+2
     j    =  blkLimits(HIGH,IAXIS) - i + 1
     SS(j) = (SS(j) - LD(j)*SS(j+1))/DD(j)
  end do
  
#else    
  do  i  = blkLimits(LOW,IAXIS)+1, blkLimits(HIGH,IAXIS)
     R     = UD(i)/DD(i-1)
     DD(i) = DD(i)-R*LD(i-1)
     SS(i) = SS(i)-R*SS(i-1)
  end do
  
  ! Back Substitution.
  SS (blkLimits(HIGH,IAXIS)) = SS(blkLimits(HIGH,IAXIS))/DD(blkLimits(HIGH,IAXIS))   
  
  do i    =  2, (blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS))+1
     j    =  blkLimits(HIGH,IAXIS) - i + 1
     SS(j) = (SS(j) - LD(j)*SS(j+1))/DD(j)
  end do
#endif



  
  
end subroutine Thomas3D



SUBROUTINE IncompleteCholesky(A,N)  
  ! Subroutine overwrites the sparse SPD matrix A
  ! with L from its IC factorization LL^T
  
  Implicit None
  
  integer, intent (IN)                    :: N
  real,    intent (INOUT), Dimension(N,N) :: A

  integer :: i,j,k

!!$  print*, "before IC"
!!$  do i=1, N
!!$     write (*,'(1p8e14.6)') A(i,1:N)
!!$  end do  
  
  do k=1,N
     A(k,k)=sqrt(A(k,k))
     do i=k+1,N
        If (A(i,k).ne.0) Then
           A(i,k)=A(i,k)/A(k,k)
        End If
     End Do
     Do j=k+1,N
        Do i=j,N
           If (A(i,j).ne.0) Then
              A(i,j)=A(i,j)-A(i,k)*A(j,k)
           End If
        End Do
     End Do
  End Do
!!$  
!!$  print*, "after IC"
!!$  do i=1, N
!!$     write (*,'(1p8e14.6)') A(i,1:N)
!!$  end do

end SUBROUTINE IncompleteCholesky





