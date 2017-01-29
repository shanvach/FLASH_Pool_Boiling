! Subroutine outtotecplot
!
! Subroutine to write out to Tecplot data in binary form.
!
! ---------------------------------------------------------------------------
#include "constants.h"
#include "Flash.h"

  subroutine sm_setStokes(blockList,blockCount,niu,w,z0)


  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkPtr, &
                 Grid_releaseBlkPtr, Grid_getBlkIndexLimits, &
                 Grid_getBlkBoundBox,Grid_getBlkCenterCoords

  use ins_interface, only : ins_velgradtensor

#ifdef FLASH_GRID_PARAMESH
  use physicaldata, only : interp_mask_unk, interp_mask_unk_res
#endif
  implicit none
#include "Flash_mpi.h"
  integer, intent(in) :: blockCount
  integer, intent(in) :: blockList(MAXBLOCKS)
  real, intent(in) :: niu, w, z0

  ! Local variables    
  integer :: numblocks,lb, blockID

  real xedge(NXB+1),xcell(NXB+1)
  real yedge(NYB+1),ycell(NYB+1)
  real zedge(NZB+1),zcell(NZB+1)
  real intsx(NXB+1), intsy(NYB+1), intsz(NZB+1)

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  real del(3),dx,dy,dz
  real, dimension(MDIM)  :: coord,bsize
  real ::  boundBox(2,MDIM)

  real :: a, r, theta, phi, x, y, z, pi
  real :: vr, vt, vx, vy, vz
  real :: xe0, ye0, ze0, xc0, yc0, zc0
 
  integer :: i, j, k

  a = 0.5
  pi = acos(-1.0)


  do lb = 1,blockcount
  
     blockID =  blockList(lb)      

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)
     dx = del(IAXIS)
     dy = del(JAXIS)
     dz = del(KAXIS)

     ! Get Coord and Bsize for the block:
     ! Bounding box:
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)

     call Grid_getBlkCenterCoords(blockId,coord)

     xe0 = coord(IAXIS) - bsize(IAXIS)/2.0 - NGUARD*dx
     xc0 = xe0 + dx/2.0

     ye0 = coord(JAXIS) - bsize(JAXIS)/2.0 - NGUARD*dy
     yc0 = ye0 + dy/2.0

     ze0 = coord(KAXIS) - bsize(KAXIS)/2.0 - NGUARD*dz - z0
     zc0 = ze0 + dz/2.0

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     do k = GRID_KLO_GC, GRID_KHI_GC
        do j = GRID_JLO_GC, GRID_JHI_GC
          do i = GRID_ILO_GC, GRID_IHI_GC
            x = xc0 + (i-1)*dx
            y = yc0 + (j-1)*dy
            z = zc0 + (k-1)*dz
            r = sqrt(x*x+y*y+z*z)
            if( abs(r) < 0.5) r = 1.0d14
            theta = acos(z/r)
            solnData(PRES_VAR,i,j,k) = -1.5*niu*w*a*cos(theta)/(r*r*r)  
          enddo
        enddo
     enddo

     do k = GRID_KLO_GC, GRID_KHI_GC
        do j = GRID_JLO_GC, GRID_JHI_GC
          do i = GRID_ILO_GC, GRID_IHI_GC+1
            x = xe0 + (i-1)*dx
            y = yc0 + (j-1)*dy
            z = zc0 + (k-1)*dz
            r = sqrt(x*x+y*y+z*z)
            if( abs(r) < 0.5) r = 1.0d14
            theta = acos(z/r)
            if(x > 1.0d-12) then
              if( y > 1.0d-12) then
                phi = atan(y/x)
              elseif(y < -1.0d-12) then
                phi = atan(y/x) + 2.0*pi
              else
                phi = 0.
              endif
            elseif( x < -1.0d-12) then
              phi = atan(y/x) + pi
            else
              if(y > 0.0) then
                phi = pi/2.0
              else
                phi = 3.0*pi/2.0
              endif
            endif
            vr = w*cos(theta)*(1+a*a*a/(2.0*r*r*r)-3.0*a/(2.0*r))
            vt =-w*sin(theta)*(1-a*a*a/(4.0*r*r*r)-3.0*a/(4.0*r))
            vx = vr*sin(theta)*cos(phi) + vt*cos(theta)*cos(phi)
            facexData(VELC_FACE_VAR,i,j,k) = vx  
          enddo
        enddo
     enddo

     do k = GRID_KLO_GC, GRID_KHI_GC
        do j = GRID_JLO_GC, GRID_JHI_GC+1
          do i = GRID_ILO_GC, GRID_IHI_GC
            x = xc0 + (i-1)*dx
            y = ye0 + (j-1)*dy
            z = zc0 + (k-1)*dz
            r = sqrt(x*x+y*y+z*z)
            if( abs(r) < 0.5) r = 1.0d14
            theta = acos(z/r)
            if(x > 1.0d-12) then
              if( y > 1.0d-12) then
                phi = atan(y/x)
              elseif(y < -1.0d-12) then
                phi = atan(y/x) + 2.0*pi
              else
                phi = 0.
              endif
            elseif( x < -1.0d-12) then
              phi = atan(y/x) + pi
            else
              if(y > 0.0) then
                phi = pi/2.0
              else
                phi = 3.0*pi/2.0
              endif
            endif
            vr = w*cos(theta)*(1+a*a*a/(2.0*r*r*r)-3.0*a/(2.0*r))
            vt =-w*sin(theta)*(1-a*a*a/(4.0*r*r*r)-3.0*a/(4.0*r))
            vy = vr*sin(theta)*sin(phi) + vt*cos(theta)*sin(phi)
            faceyData(VELC_FACE_VAR,i,j,k) = vy  
          enddo
        enddo
     enddo

     do k = GRID_KLO_GC, GRID_KHI_GC+1
        do j = GRID_JLO_GC, GRID_JHI_GC
          do i = GRID_ILO_GC, GRID_IHI_GC
            x = xc0 + (i-1)*dx
            y = yc0 + (j-1)*dy
            z = ze0 + (k-1)*dz
            r = sqrt(x*x+y*y+z*z)
            if( abs(r) < 0.5) r = 1.0d14
            theta = acos(z/r)
            if(x > 1.0d-12) then
              if( y > 1.0d-12) then
                phi = atan(y/x)
              elseif(y < -1.0d-12) then
                phi = atan(y/x) + 2.0*pi
              else
                phi = 0.
              endif
            elseif( x < -1.0d-12) then
              phi = atan(y/x) + pi
            else
              if(y > 0.0) then
                phi = pi/2.0
              else
                phi = 3.0*pi/2.0
              endif
            endif
            vr = w*cos(theta)*(1+a*a*a/(2.0*r*r*r)-3.0*a/(2.0*r))
            vt =-w*sin(theta)*(1-a*a*a/(4.0*r*r*r)-3.0*a/(4.0*r))
            vz = vr*cos(theta) - vt*sin(theta)
            facezData(VELC_FACE_VAR,i,j,k) = vz  
          enddo
        enddo
     enddo

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  enddo

  return
  End subroutine sm_setStokes

