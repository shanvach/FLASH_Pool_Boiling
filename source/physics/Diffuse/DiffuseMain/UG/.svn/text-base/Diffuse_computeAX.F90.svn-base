!!****if* source/physics/Diffuse/DiffuseMain/UG/Diffuse_computeAX
!!
!!  NAME 
!!
!!  diff_computeAX
!!
!!  SYNOPSIS
!!
!!  diff_computeAX (integer, intent(IN)           :: blkLimits
!!                  integer, intent(IN)           :: blkLimitsGC
!!                  real,    intent(IN)           :: X    
!!                  real,    intent(OUT)          :: AX 
!!                  real,    intent(IN)           :: iFactorB  
!!                  real,    intent(IN)           :: iFactorA    
!!                  real,    intent(IN), OPTIONAL :: iFactorC 
!!                  real,    intent(IN)           :: del
!!                  real,    intent(IN)           :: dt
!!                  real,    intent(IN)           :: theta)
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
!! NOTES
!!  
!!
!!***

subroutine Diffuse_computeAX (X, AX, iFactorB, iFactorA, blkLimits, blkLimitsGC,del, dt, theta, centerCoords,iFactorC)
  
  use Diffuse_data,   ONLY : diff_geometry

  implicit none 

#include "Flash.h"
#include "constants.h"

!!------------------------------------------------------------------------------------------
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC  
  real, intent(IN)          :: X  (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))  
  real, intent(OUT)         :: AX (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
  real, intent(IN)          :: iFactorB (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) 
  real, intent(IN)          :: iFactorA (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))  
  real,intent(IN), OPTIONAL :: iFactorC (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
  real,intent(IN)           :: centerCoords (blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1)
  
  
  real, dimension(MDIM), intent(IN) :: del
  real, intent (IN)                 :: dt
  real, intent (IN)                 :: theta 
  
  !!-------------------------------------------------------------------------------------------
  integer               :: i,j,k  
  real                  :: CondR, CondL
  real                  :: DD, LD, UD, CDiv
  !! CYLINDRICAL
  real                  :: r, dr, rphalf, rmhalf, dxsq
  
  AX = 0.0    

  if (diff_geometry == CARTESIAN) then 
     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
              
              CDiv  = 1.0/iFactorA(i,j,k)
              
              CondR = 0.5 *  (iFactorB(i+1,j,k)*CDiv + iFactorB(i,j,k)*CDiv)
              CondL = 0.5 *  (iFactorB(i-1,j,k)*CDiv + iFactorB(i,j,k)*CDiv)
              
              LD = -theta*CondL*(dt/(del(1)**2))
              UD = -theta*CondR*(dt/(del(1)**2))
              DD = 1.0 + theta*(CondR + CondL)*(dt/(del(1)**2))
              
              AX(i,j,k) = LD*X(i-1,j,k) + UD*X(i+1,j,k)           
              
#if NDIM >= 2
              CondR = 0.5 *  CDiv * (iFactorB(i,j+1,k) + iFactorB(i,j,k))
              CondL = 0.5 *  CDiv * (iFactorB(i,j-1,k) + iFactorB(i,j,k))
              
              LD = -theta*CondL*(dt/(del(2)**2))
              UD = -theta*CondR*(dt/(del(2)**2))
              
              AX(i,j,k) = AX(i,j,k) + (LD*X(i,j-1,k) + UD*X(i,j+1,k))
              
              DD = DD + theta*(CondR + CondL)*(dt/(del(2)**2))                      
#if NDIM == 3
           
              CondR = 0.5 *  CDiv * (iFactorB(i,j,k+1) + iFactorB(i,j,k))
              CondL = 0.5 *  CDiv * (iFactorB(i,j,k-1) + iFactorB(i,j,k))
              
              LD = -theta*CondL*(dt/(del(3)**2))
              UD = -theta*CondR*(dt/(del(3)**2))
              
              AX(i,j,k) = AX(i,j,k) + (LD*X(i,j,k-1) + UD*X(i,j,k+1))          
           
              DD = DD + theta*(CondR + CondL)*(dt/(del(3)**2))
#endif          
              
#endif                  
              if (present(iFactorC)) then
                 DD = DD + dt * iFactorC (i,j,k)
              end if
              
              AX(i,j,k) = AX(i,j,k) + DD*X(i,j,k)
              
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
              
              CDiv  = 1.0/iFactorA(i,j,k)
              
              CondR = 0.5 *  (iFactorB(i+1,j,k)*CDiv + iFactorB(i,j,k)*CDiv)
              CondL = 0.5 *  (iFactorB(i-1,j,k)*CDiv + iFactorB(i,j,k)*CDiv)           
           
              !! The matrix is modified.
              LD = -theta*CondL*rmhalf*dt/dxsq
              UD = -theta*CondR*rphalf*dt/dxsq
              DD = 1.0 + theta*(CondR*rphalf + CondL*rmhalf)*dt/dxsq
              
              AX(i,j,k) = LD*X(i-1,j,k) + UD*X(i+1,j,k)              
           
#if NDIM >= 2              
              !! in Cylindrical this would be Z-direction and is handled the same way as before. 
              
              CondR = 0.5 *  CDiv * (iFactorB(i,j+1,k) + iFactorB(i,j,k))
              CondL = 0.5 *  CDiv * (iFactorB(i,j-1,k) + iFactorB(i,j,k))
              
              LD = -theta*CondL*(dt/(del(2)**2))
              UD = -theta*CondR*(dt/(del(2)**2))
           
              AX(i,j,k) = AX(i,j,k) + (LD*X(i,j-1,k) + UD*X(i,j+1,k))          
              
              DD = DD + theta*(CondR + CondL)*(dt/(del(3)**2))
              
#if NDIM == 3
              call Driver_abortFlash("CYLINDRICAL geometry : PHI direction not yet supported!")              
#endif          
              
#endif                
              if (present(iFactorC)) then
                 DD = DD + dt * iFactorC (i,j,k)
              end if
              
              AX(i,j,k) = AX(i,j,k) + DD*X(i,j,k)
              
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
              
              CDiv  = 1.0/iFactorA(i,j,k)
              
              CondR = 0.5 *  (iFactorB(i+1,j,k)*CDiv + iFactorB(i,j,k)*CDiv)
              CondL = 0.5 *  (iFactorB(i-1,j,k)*CDiv + iFactorB(i,j,k)*CDiv)           
              
              !! The matrix is modified.
              LD = -theta*CondL*(rmhalf**2)*dt/dxsq
              UD = -theta*CondR*(rphalf**2)*dt/dxsq
              DD = 1.0 + theta*(CondR*(rphalf**2) + CondL*(rmhalf**2))*dt/dxsq
              
              AX(i,j,k) = LD*X(i-1,j,k) + UD*X(i+1,j,k)
              
              if (present(iFactorC)) then
                 DD = DD + dt * iFactorC (i,j,k)
              end if
              
              AX(i,j,k) = AX(i,j,k) + DD*X(i,j,k)
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
 
end subroutine Diffuse_computeAX

