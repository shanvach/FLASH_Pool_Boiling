!!****if* source/physics/Diffuse/DiffuseMain/UG/diff_computeB
!!
!!  NAME 
!!
!!  diff_computeILU
!!
!!  SYNOPSIS
!!
!!  call diff_computeB
!!
!!
!!  DESCRIPTION 
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

subroutine diff_computeB (blockCount, blockList, iVar, iFactorA, iFactorB, dt, theta, bcTypes, bcValues, iFactorD)
  
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Diffuse_data,     ONLY : diff_geometry, diff_speedlt
  use Grid_interface,   ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,    &
                               Grid_getBlkIndexLimits,Grid_getDeltas, &
                               Grid_getCellCoords, &
                               Grid_getBlkBC, &
                               GRID_PDE_BND_PERIODIC,  &
                               GRID_PDE_BND_NEUMANN,   &
                               GRID_PDE_BND_DIRICHLET
  implicit none
  
#include "Flash.h"
#include "constants.h" 
  
  integer, intent(IN) :: iVar
  integer, intent(IN) :: iFactorB
  integer, intent(IN) :: iFactorA
  integer, OPTIONAL, intent(IN) :: iFactorD
  real, intent(IN) :: dt
  real, intent(IN) :: theta
  integer, intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  integer, intent(IN) :: bcTypes(6)
  real,    intent(IN) :: bcValues(2,6)
  
  integer :: i, j, k, lb
  real, dimension(MDIM) :: del    
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  real :: CDiv 

  !! Needed if Geometry is cylindrical, spherical.
  real,allocatable,dimension(:) :: centerCoords
  logical :: gcell = .true. 
  real :: rad, dr, rphalf, rmhalf, dxsq, Cond_R, Cond_L, c

  integer, dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC  
  integer :: datasize
  integer:: faces(2,MDIM)
  
  call Timers_start("diff_computeB")  

  !! To store center of cell, needed only for Spherical, Cylindrical.
  if (diff_geometry /= CARTESIAN) then
     call Grid_getBlkIndexLimits(blockList(blockCount),blkLimits,blkLimitsGC) 
     datasize = blkLimitsGC(HIGH, IAXIS)- blkLimitsGC(LOW, IAXIS)+1
     allocate(centerCoords(datasize))
  end if
  
  
  do lb = 1, blockCount    
     call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC) 
     call Grid_getDeltas(blocklist(lb), del)    
     call Grid_getBlkPtr(blocklist(lb), solnVec)  
     call Grid_getBlkBC (blockList(lb), faces)
     

     !! APPLY BC - Should move to Grid_applyBC 
     if (faces(1,IAXIS) /= NOT_BOUNDARY) then
        if (bcTypes(1) == GRID_PDE_BND_NEUMANN) then
           solnVec(iFactorB,blkLimits(LOW,IAXIS)-1,:,:)  = solnVec(iFactorB,blkLimits(LOW,IAXIS),:,:)
           solnVec(iVar,    blkLimits(LOW,IAXIS)-1,:,:)  = solnVec(iVar    ,blkLimits(LOW,IAXIS),:,:)
        else if (bcTypes(1) == VACUUM) then         
           c = 2.0/(diff_speedlt*del(1))
           solnVec(iFactorB,blkLimits(LOW,IAXIS)-1,:,:)  = solnVec(iFactorB,blkLimits(LOW,IAXIS),:,:)
           solnVec(iVar,    blkLimits(LOW,IAXIS)-1,:,:)  = solnVec(iVar    ,blkLimits(LOW,IAXIS),:,:)* &
               (c*solnVec(iFactorB,blkLimits(LOW,IAXIS),:,:) - 0.5)/(c*solnVec(iFactorB,blkLimits(LOW,IAXIS),:,:) + 0.5)
        else if (bcTypes(1) == GRID_PDE_BND_DIRICHLET) then
           solnVec(iFactorB,blkLimits(LOW,IAXIS)-1,:,:) =  2.0*bcValues(2,1) - solnVec(iFactorB,blkLimits(LOW,IAXIS),:,:)!! 
           solnVec(iVar,    blkLimits(LOW,IAXIS)-1,:,:) =  (2.0*bcValues(1,1) - solnVec(iVar,    blkLimits(LOW,IAXIS),:,:))
        end if
     end if

     if (faces(2,IAXIS) /= NOT_BOUNDARY) then
        if (bcTypes(2) == GRID_PDE_BND_NEUMANN) then
           solnVec(iFactorB,blkLimits(HIGH,IAXIS)+1,:,:)  = solnVec(iFactorB,blkLimits(HIGH,IAXIS),:,:)
           solnVec(iVar,    blkLimits(HIGH,IAXIS)+1,:,:)  = solnVec(iVar    ,blkLimits(HIGH,IAXIS),:,:)
        else if (bcTypes(2) == VACUUM) then           
           c = 2.0/(diff_speedlt*del(1))
           solnVec(iFactorB,blkLimits(HIGH,IAXIS)+1,:,:)  = solnVec(iFactorB,blkLimits(HIGH,IAXIS),:,:)
           solnVec(iVar,    blkLimits(HIGH,IAXIS)+1,:,:)  = solnVec(iVar    ,blkLimits(HIGH,IAXIS),:,:)* &
               (c*solnVec(iFactorB,blkLimits(HIGH,IAXIS),:,:) + 0.5)/(c*solnVec(iFactorB,blkLimits(HIGH,IAXIS),:,:) - 0.5)
        else if (bcTypes(2) == GRID_PDE_BND_DIRICHLET) then
           solnVec(iFactorB,blkLimits(HIGH,IAXIS)+1,:,:) =  2.0*bcValues(2,2) - solnVec(iFactorB,blkLimits(HIGH,IAXIS),:,:)
           solnVec(iVar,    blkLimits(HIGH,IAXIS)+1,:,:) =  2.0*bcValues(1,2) - solnVec(iVar,    blkLimits(HIGH,IAXIS),:,:)
        end if
     end if
     
#if NDIM >= 2
     if (faces(1,JAXIS) /= NOT_BOUNDARY) then
        if (bcTypes(3) == GRID_PDE_BND_NEUMANN) then
           solnVec(iFactorB,:,blkLimits(LOW,JAXIS)-1,:)  = solnVec(iFactorB,:,blkLimits(LOW,JAXIS),:)
           solnVec(iVar,    :,blkLimits(LOW,JAXIS)-1,:)  = solnVec(iVar    ,:,blkLimits(LOW,JAXIS),:)
        else if (bcTypes(3) == VACUUM) then
           c = 2.0/(diff_speedlt*del(2))
           solnVec(iFactorB,:,blkLimits(LOW,JAXIS)-1,:)  = solnVec(iFactorB,:,blkLimits(LOW,JAXIS),:)
           solnVec(iVar,    :,blkLimits(LOW,JAXIS)-1,:)  = solnVec(iVar    ,:,blkLimits(LOW,JAXIS),:)* &
               (c*solnVec(iFactorB,:,blkLimits(LOW,JAXIS),:) - 0.5)/(c*solnVec(iFactorB,:,blkLimits(LOW,JAXIS),:) + 0.5)
        else if (bcTypes(3) == GRID_PDE_BND_DIRICHLET) then
           solnVec(iFactorB,:,blkLimits(LOW,JAXIS)-1,:) =  2.0*bcValues(2,3) - solnVec(iFactorB,:,blkLimits(LOW,JAXIS),:)
           solnVec(iVar,    :,blkLimits(LOW,JAXIS)-1,:) =  2.0*bcValues(1,3) - solnVec(iVar,    :,blkLimits(LOW,JAXIS),:)
        end if
     end if
     
     if (faces(2,JAXIS) /= NOT_BOUNDARY) then
        if (bcTypes(4) == GRID_PDE_BND_NEUMANN) then
           solnVec(iFactorB,:,blkLimits(HIGH,JAXIS)+1,:)  = solnVec(iFactorB,:,blkLimits(HIGH,JAXIS),:)
           solnVec(iVar,    :,blkLimits(HIGH,JAXIS)+1,:)  = solnVec(iVar    ,:,blkLimits(HIGH,JAXIS),:)
        else if (bcTypes(4) == VACUUM) then
           c = 2.0/(diff_speedlt*del(2))
           solnVec(iFactorB,:,blkLimits(HIGH,JAXIS)+1,:)  = solnVec(iFactorB,:,blkLimits(HIGH,JAXIS),:)
           solnVec(iVar,    :,blkLimits(HIGH,JAXIS)+1,:)  = solnVec(iVar    ,:,blkLimits(HIGH,JAXIS),:)* &
               (c*solnVec(iFactorB,:,blkLimits(HIGH,JAXIS),:) + 0.5)/(c*solnVec(iFactorB,:,blkLimits(HIGH,JAXIS),:) - 0.5)
        else if (bcTypes(4) == GRID_PDE_BND_DIRICHLET) then
           solnVec(iFactorB,:,blkLimits(HIGH,JAXIS)+1,:) =  2.0*bcValues(2,4) - solnVec(iFactorB,:,blkLimits(HIGH,JAXIS),:)
           solnVec(iVar,    :,blkLimits(HIGH,JAXIS)+1,:) =  2.0*bcValues(1,4) - solnVec(iVar,    :,blkLimits(HIGH,JAXIS),:)           
        end if
     end if
     
     
#if NDIM == 3
     if (faces(1,KAXIS) /= NOT_BOUNDARY) then
        if (bcTypes(5) == GRID_PDE_BND_NEUMANN) then
           solnVec(iFactorB,:,:,blkLimits(LOW,KAXIS)-1)  = solnVec(iFactorB,:,:,blkLimits(LOW,KAXIS))
           solnVec(iVar,    :,:,blkLimits(LOW,KAXIS)-1)  = solnVec(iVar    ,:,:,blkLimits(LOW,KAXIS))
        else if (bcTypes(5) == VACUUM) then
           c = 2.0/(diff_speedlt*del(3))
           solnVec(iFactorB,:,:,blkLimits(LOW,KAXIS)-1)  = solnVec(iFactorB,:,:,blkLimits(LOW,KAXIS))
           solnVec(iVar,    :,:,blkLimits(LOW,KAXIS)-1)  = solnVec(iVar    ,:,:,blkLimits(LOW,KAXIS))* &
               (c*solnVec(iFactorB,:,:,blkLimits(LOW,KAXIS)) - 0.5)/(c*solnVec(iFactorB,:,:,blkLimits(LOW,KAXIS)) + 0.5)
        else if (bcTypes(5) == GRID_PDE_BND_DIRICHLET) then
           solnVec(iFactorB,:,:,blkLimits(LOW,KAXIS)-1) =  2.0*bcValues(2,5) - solnVec(iFactorB,:,:,blkLimits(LOW,KAXIS))
           solnVec(iVar,    :,:,blkLimits(LOW,KAXIS)-1) =  2.0*bcValues(1,5) - solnVec(iVar,    :,:,blkLimits(LOW,KAXIS))          
        end if
     end if
     
     if (faces(2,KAXIS) /= NOT_BOUNDARY) then
        if (bcTypes(6) == GRID_PDE_BND_NEUMANN) then
           solnVec(iFactorB,:,:,blkLimits(HIGH,KAXIS)+1)  = solnVec(iFactorB,:,:,blkLimits(HIGH,KAXIS))
           solnVec(iVar,    :,:,blkLimits(HIGH,KAXIS)+1)  = solnVec(iVar    ,:,:,blkLimits(HIGH,KAXIS))
        else if (bcTypes(6) == VACUUM) then
           c = 2.0/(diff_speedlt*del(3))
           solnVec(iFactorB,:,:,blkLimits(HIGH,KAXIS)+1)  = solnVec(iFactorB,:,:,blkLimits(HIGH,KAXIS))
           solnVec(iVar,    :,:,blkLimits(HIGH,KAXIS)+1)  = solnVec(iVar    ,:,:,blkLimits(HIGH,KAXIS),:)* &
               (c*solnVec(iFactorB,:,:,blkLimits(HIGH,KAXIS)) + 0.5)/(c*solnVec(iFactorB,:,:,blkLimits(HIGH,KAXIS)) - 0.5)
        else if (bcTypes(6) == GRID_PDE_BND_DIRICHLET) then
           solnVec(iFactorB,:,:,blkLimits(HIGH,KAXIS)+1) =  2.0*bcValues(2,6)-solnVec(iFactorB,:,:,blkLimits(HIGH,KAXIS))
           solnVec(iVar,    :,:,blkLimits(HIGH,KAXIS)+1) =  2.0*bcValues(1,6)-solnVec(iVar,    :,:,blkLimits(HIGH,KAXIS))
           
        end if
     end if
#endif
     

#endif    
     

     
     !!===============================================================================================                     
     !! CARTESIAN 1D, 2D, 3D.
     !!===============================================================================================                     
     if (diff_geometry == CARTESIAN) then 
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 cDiv = 1.0!/solnVec(iFactorA,i,j,k)
                 
                 Cond_R = 0.5*(solnVec(iFactorB,i+1,j,k) + solnVec(iFactorB,i,j,k))*CDiv
                 Cond_L = 0.5*(solnVec(iFactorB,i-1,j,k) + solnVec(iFactorB,i,j,k))*CDiv
                 
                 SolnVec(WTMP_VAR,i,j,k) = solnVec(iFactorA,i,j,k)*SolnVec(iVar,i,j,k) + (1.0-theta)*(dt/del(1)**2) * &
                      (cond_R*solnVec(iVar,i+1,j,k) - (cond_R+cond_L)*solnVec(iVar,i,j,k) + cond_L*solnVec(iVar,i-1,j,k))
                 
#if NDIM >= 2
                 Cond_R = 0.5*(solnVec(iFactorB,i,j+1,k) + solnVec(iFactorB,i,j,k))*CDiv
                 Cond_L = 0.5*(solnVec(iFactorB,i,j-1,k) + solnVec(iFactorB,i,j,k))*CDiv
                 
                 SolnVec(WTMP_VAR,i,j,k) =  SolnVec(WTMP_VAR,i,j,k) + (1.0-theta)*(dt/del(2)**2) * &
                      (cond_R*solnVec(iVar,i,j+1,k)-(cond_R+cond_L)*solnVec(iVar,i,j,k)+cond_L*solnVec(iVar,i,j-1,k))



#if NDIM == 3
                 Cond_R = 0.5*(solnVec(iFactorB,i,j,k+1) + solnVec(iFactorB,i,j,k))*CDiv
                 Cond_L = 0.5*(solnVec(iFactorB,i,j,k-1) + solnVec(iFactorB,i,j,k))*CDiv
                 
                 SolnVec(WTMP_VAR,i,j,k) =  (SolnVec(WTMP_VAR,i,j,k) + (1.0-theta)*(dt/del(3)**2) * &
                      (cond_R*solnVec(iVar,i,j,k+1) - (cond_R+cond_L)*solnVec(iVar,i,j,k) + cond_L*solnVec(iVar,i,k-1,k)))
#endif              
#endif              

                 if (present(iFactorD)) then
                    SolnVec(WTMP_VAR,i,j,k) = SolnVec(WTMP_VAR,i,j,k) + dt*SolnVec(iFactorD, i,j,k)
                 end if
                 
                 
             end do
          end do
       end do       
       
       
       !!===============================================================================================                            
       !! Cylindrical 1D, 2D.
       !!===============================================================================================                     
    else if (diff_geometry == CYLINDRICAL) then
       call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, .true., centerCoords, blkLimits(HIGH, IAXIS))

       do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
          do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
             do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                rad  = centerCoords(i)           
                dr = del(1) 
                rphalf = rad + 0.5*dr
                rmhalf = rad - 0.5*dr           
                dxsq = rad*dr**2
                
                cDiv = 1.0!/solnVec(iFactorA,i,j,k)

                
                Cond_R = 0.5*(solnVec(iFactorB,i+1,j,k) + solnVec(iFactorB,i,j,k))*CDiv
                Cond_L = 0.5*(solnVec(iFactorB,i-1,j,k) + solnVec(iFactorB,i,j,k))*CDiv
                
                SolnVec(WTMP_VAR,i,j,k) = &
                     solnVec(iFactorA,i,j,k)*SolnVec(iVar,i,j,k) + (1.0-theta)*(dt/dxsq)*(Cond_L*(rmhalf)*solnVec(iVar,i-1,j,k) & 
                     -(Cond_L*(rmhalf)+Cond_R*(rphalf))*solnVec(iVar,i,j,k) &
                     + Cond_R*(rphalf)*solnVec(iVar,i+1,j,k))




#if NDIM >= 2
                Cond_R = 0.5*(solnVec(iFactorB,i,j+1,k) + solnVec(iFactorB,i,j,k))*CDiv
                Cond_L = 0.5*(solnVec(iFactorB,i,j-1,k) + solnVec(iFactorB,i,j,k))*CDiv
                
                SolnVec(WTMP_VAR,i,j,k) =  SolnVec(WTMP_VAR,i,j,k) + ((1.0-theta)*(dt/del(2)**2) * &
                     (cond_R*solnVec(iVar,i,j+1,k)-(cond_R+cond_L)*solnVec(iVar,i,j,k)+cond_L*solnVec(iVar,i,j-1,k)))
#endif
                
                if (present(iFactorD)) then
                   SolnVec(WTMP_VAR,i,j,k) = SolnVec(WTMP_VAR,i,j,k) + dt*SolnVec(iFactorD, i,j,k)
                end if
                
                SolnVec(WTMP_VAR,i,j,k) =  rad*dr*SolnVec(WTMP_VAR,i,j,k) !! in R-Z RHS is scaled by r*dr to keep A symmetric                
                
                
             end do
          end do
       end do
       
       !!===============================================================================================                     
       !!Spherical 1D.
       !!===============================================================================================                     
    else if (diff_geometry == SPHERICAL) then            
       call Grid_getCellCoords(IAXIS, lb, CENTER, gcell,centerCoords, datasize)
       do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
          do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
             do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                
                rad  = centerCoords(i)           
                dr = del(1)           
             
                rphalf = rad + 0.5*dr
                rmhalf = rad - 0.5*dr           
                dxsq = (rad**2)*(dr**2)
             
                cDiv = 1.0/solnVec(iFactorA,i,j,k)
                    
                Cond_R = 0.5*(solnVec(iFactorB,i+1,j,k) + solnVec(iFactorB,i,j,k))*CDiv
                Cond_L = 0.5*(solnVec(iFactorB,i-1,j,k) + solnVec(iFactorB,i,j,k))*CDiv
            
                SolnVec(WTMP_VAR,i,j,k) = &
                     SolnVec(iVar,i,j,k) + (1.0-theta)*(dt/dxsq)*(Cond_L*(rmhalf**2)*solnVec(iVar,i-1,j,k) & 
                                                                -(Cond_L*(rmhalf**2)+Cond_R*(rphalf**2))*solnVec(iVar,i,j,k) &
                                                                + Cond_R*(rphalf**2)*solnVec(iVar,i+1,j,k))           
             end do
          end do
       end do
    end if
    

    
    
    call Grid_releaseBlkPtr(blocklist(lb), solnVec)  
 end do
 !!===============================================================================================

 if (diff_geometry /= CARTESIAN) then
    deallocate (centerCoords)
 end if


 call Timers_stop("diff_computeB")   

end subroutine diff_computeB
