!!****if* source/physics/Diffuse/DiffuseMain/UG/diff_computeAX
!!
!!  NAME 
!!
!!  diff_computeAX
!!
!!  SYNOPSIS
!!
!!  call diff_computeAX (integer(IN)           :: blockID,
!!                  integer(IN)           :: blkLimits,
!!                  integer(IN)           :: blkLimitsGC,
!!                  real(IN)              :: X    ,
!!                  real(OUT)             :: AX ,
!!                  real(IN)              :: iFactorA    ,
!!                  real(IN)              :: iFactorB ,
!!                  real(IN)   , OPTIONAL :: iFactorC ,
!!                  real(IN)              :: dt,
!!                  real(IN)              :: theta)
!!
!!
!!  DESCRIPTION 
!!      This routine computes AX, Jacobian free, as we do not have to store AX.
!!      The elements of A are built from FD of general diffusion/conduction equation      
!!
!!      A * dV/dt = d/dx(B*dV/dx) + d/dy(B*dV/dy) + d/dx(B*dV/dz) + C*V + D
!!
!!
!! ARGUMENTS
!!      
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
subroutine diff_computeAX (blockID, X, AX, iFactorB, iFactorA, blkLimits, blkLimitsGC, dt, theta,iFactorC)
           
  use Diffuse_data,   ONLY : diff_geometry
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &                                
       Grid_getDeltas
  use Timers_interface, ONLY : Timers_start, Timers_stop
  
  implicit none 
  
#include "Flash.h"
#include "constants.h"
  
  !!------------------------------------------------------------------------------------------
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC  
  integer,intent(IN) :: X 
  integer,intent(IN) :: iFactorB 
  integer,intent(IN) :: iFactorA  
  real,intent(OUT) :: AX (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                          blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                          blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) 
  integer, intent(IN):: blockID
  real, intent (IN)  :: dt
  real, intent (IN)  :: theta   
  integer,intent(IN),OPTIONAL :: iFactorC

  
  !!-------------------------------------------------------------------------------------------
  integer               :: i,j,k  
  real                  :: CondR, CondL
  real                  :: DD, LD, UD, CDiv
  !! CYLINDRICAL
  real                  :: r, dr, rphalf, rmhalf, dxsq

  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  real, dimension(MDIM) :: del
  real, allocatable, dimension(:) :: centerCoords
  integer :: datasize


  call Timers_start("diff_computeAX")    

  !! To store center of cell, needed only for Spherical, Cylindrical.
  if (diff_geometry /= CARTESIAN) then
     datasize = blkLimitsGC(HIGH, IAXIS)- blkLimitsGC(LOW, IAXIS)+1
     allocate(centerCoords(datasize))
  end if
  
  call Grid_getDeltas(blockID, del)
  call Grid_getBlkPtr(blockID, solnVec)   
  
  AX = 0.0    
  
  if (diff_geometry == CARTESIAN) then 
     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
              
              CDiv  = 1.0/solnVec(iFactorA,i,j,k)
              
              CondR = 0.5 *  (solnVec(iFactorB,i+1,j,k)*CDiv + solnVec(iFactorB,i,j,k)*CDiv)
              CondL = 0.5 *  (solnVec(iFactorB,i-1,j,k)*CDiv + solnVec(iFactorB,i,j,k)*CDiv)
              
              LD = -theta*CondL*(dt/(del(1)**2))
              UD = -theta*CondR*(dt/(del(1)**2))
              DD = 1.0 + theta*(CondR + CondL)*(dt/(del(1)**2))
              
              AX(i,j,k) = LD*solnVec(X,i-1,j,k) + UD*solnVec(X,i+1,j,k)           
              
#if NDIM >= 2
              CondR = 0.5 *  CDiv * (solnVec(iFactorB,i,j+1,k) + solnVec(iFactorB,i,j,k))
              CondL = 0.5 *  CDiv * (solnVec(iFactorB,i,j-1,k) + solnVec(iFactorB,i,j,k))
              
              LD = -theta*CondL*(dt/(del(2)**2))
              UD = -theta*CondR*(dt/(del(2)**2))
              
              AX(i,j,k) = AX(i,j,k) + (LD*solnVec(X,i,j-1,k) + UD*solnVec(X,i,j+1,k))
              
              DD = DD + theta*(CondR + CondL)*(dt/(del(2)**2))                      
#if NDIM == 3
           
              CondR = 0.5 *  CDiv * (solnVec(iFactorB,i,j,k+1) + solnVec(iFactorB,i,j,k))
              CondL = 0.5 *  CDiv * (solnVec(iFactorB,i,j,k-1) + solnVec(iFactorB,i,j,k))
              
              LD = -theta*CondL*(dt/(del(3)**2))
              UD = -theta*CondR*(dt/(del(3)**2))
              
              AX(i,j,k) = AX(i,j,k) + (LD*solnVec(X,i,j,k-1) + UD*solnVec(X,i,j,k+1))          
           
              DD = DD + theta*(CondR + CondL)*(dt/(del(3)**2))
#endif          
              
#endif                  
              if (present(iFactorC)) then
                 DD = DD + dt * solnVec(iFactorC,i,j,k)
              end if
              
              AX(i,j,k) = AX(i,j,k) + DD*solnVec(X,i,j,k)
              
          end do
        end do
     end do
     
  else if (diff_geometry == CYLINDRICAL) then   
     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)                              
              r  = centerCoords(i)           
              dr = del(1)           
              rphalf = r + 0.5*dr
              rmhalf = r - 0.5*dr           
              dxsq = r*dr**2
              
              CDiv  = 1.0/solnVec(iFactorA,i,j,k)
              
              CondR = 0.5 *  (solnVec(iFactorB,i+1,j,k)*CDiv + solnVec(iFactorB,i,j,k)*CDiv)
              CondL = 0.5 *  (solnVec(iFactorB,i-1,j,k)*CDiv + solnVec(iFactorB,i,j,k)*CDiv)           
           
              !! The matrix is modified.
              LD = -theta*CondL*rmhalf*dt/dxsq
              UD = -theta*CondR*rphalf*dt/dxsq
              DD = 1.0 + theta*(CondR*rphalf + CondL*rmhalf)*dt/dxsq
              
              AX(i,j,k) = LD*solnVec(X,i-1,j,k) + UD*solnVec(X,i+1,j,k)              
           
#if NDIM >= 2              
              !! in Cylindrical this would be Z-direction and is handled the same way as before. 
              
              CondR = 0.5 *  CDiv * (solnVec(iFactorB,i,j+1,k) + solnVec(iFactorB,i,j,k))
              CondL = 0.5 *  CDiv * (solnVec(iFactorB,i,j-1,k) + solnVec(iFactorB,i,j,k))
              
              LD = -theta*CondL*(dt/(del(2)**2))
              UD = -theta*CondR*(dt/(del(2)**2))
           
              AX(i,j,k) = AX(i,j,k) + (LD*solnVec(X,i,j-1,k) + UD*solnVec(X,i,j+1,k))          
              
              DD = DD + theta*(CondR + CondL)*(dt/(del(3)**2))
              
#if NDIM == 3
              call Driver_abortFlash("CYLINDRICAL geometry : PHI direction not yet supported!")              
#endif          
              
#endif                
              if (present(iFactorC)) then
                 DD = DD + dt * solnVec(iFactorC,i,j,k)
              end if
              
              AX(i,j,k) = AX(i,j,k) + DD*solnVec(X,i,j,k)
              
           end do
        end do
     end do
  else if (diff_geometry == SPHERICAL) then    
     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)                             
              
#if NDIM >= 2 
              call Driver_abortFlash("SPHERICAL geometry : Only 1D supported!")
#else              
              r  = centerCoords(i)  
              
              dr = del(1)           
              
              rphalf = r + 0.5*dr
              rmhalf = r - 0.5*dr           
              dxsq = (r*dr)**2
              
              CDiv  = 1.0/solnVec(iFactorA,i,j,k)
              
              CondR = 0.5 *  (solnVec(iFactorB,i+1,j,k)*CDiv + solnVec(iFactorB,i,j,k)*CDiv)
              CondL = 0.5 *  (solnVec(iFactorB,i-1,j,k)*CDiv + solnVec(iFactorB,i,j,k)*CDiv)           
              
              !! The matrix is modified.
              LD = -theta*CondL*(rmhalf**2)*dt/dxsq
              UD = -theta*CondR*(rphalf**2)*dt/dxsq
              DD = 1.0 + theta*(CondR*(rphalf**2) + CondL*(rmhalf**2))*dt/dxsq
              
              AX(i,j,k) = LD*solnVec(X,i-1,j,k) + UD*solnVec(X,i+1,j,k)
              
              if (present(iFactorC)) then
                 DD = DD + dt * solnVec(iFactorC,i,j,k)
              end if
              
              AX(i,j,k) = AX(i,j,k) + DD*solnVec(X,i,j,k)
#endif
              
           end do
        end do
     end do
             
  else if(diff_geometry == POLAR) then
     !! Not yet implemented.
     call Driver_abortFlash("POLAR not yet supported !!")
  else
     call Driver_abortFlash("3D Cartesian, 2D Cylindrical, 1D Spherical only !")
  end if


  call Grid_releaseBlkPtr(blockID, solnVec)    
  
 if (diff_geometry /= CARTESIAN) then
    deallocate (centerCoords)
 end if

  call Timers_stop("diff_computeAX")       
  
end subroutine diff_computeAX
