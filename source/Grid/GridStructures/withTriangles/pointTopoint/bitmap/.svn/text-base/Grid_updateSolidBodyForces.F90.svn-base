!!****if* source/Grid/GridStructures/withTriangles/pointTopoint/bitmap/Grid_updateSolidBodyForces
!!
!! NAME
!!  Grid_updateSolidBodyForces
!!
!! SYNOPSIS
!!
!!  Grid_updateSolidBodyForces(integer, INTENT(in)    :: blkID,
!!                        integer, INTENT(in)    :: b,
!!                        integer, INTENT(in)    :: localParticleCount,
!!                           real, dimension(MDIM), INTENT(in)    :: particleposn)
!!  
!! DESCRIPTION 
!!  
!!  The velocity update and forcing routine
!!
!!  Overview of the algoritm
!!
!! ARGUMENTS
!!
!!  blkID:  the local block ID
!!
!!  b : solid body number
!!
!!  localParticleCount : particle number in each processor 
!!
!!  particleposn : the position coordinates of the particle
!! 
!!  
!!***

#include "constants.h"
#include "Flash.h"

!subroutine Grid_updateSolidBodyForces(blkID, b, localParticleCount, particleposn)
subroutine Grid_updateSolidBodyForces(blkID, particleData)
  use Grid_data, ONLY : gr_meshMe
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, gr_sbPtNumX, &
       gr_sbPtNumY, gr_sbPtNumZ
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkBoundBox, &
       Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr
  implicit none
  include "Flash_mpi.h"

  integer, intent(IN) :: blkID
  real, dimension(NPART_PROPS), intent(INOUT) :: particleData
!  integer, intent(IN) :: b, localParticleCount
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
!  real, dimension(MDIM), intent(IN) :: particleposn
  real, dimension(MDIM) :: particleposn
  real, dimension(MDIM) :: cellval, del
  real, allocatable, dimension(:) :: coord, dist
  integer :: coordSize, axis
  integer, dimension(MDIM) :: cell, faces, facecent
  real, pointer, dimension(:,:,:,:) :: fcx, fcy, fcz
  real :: ei, ej, ek, xcell, ycell, zcell


  ! Start of interpolation routine, from eulerian FACEVAR to particle.
        
  ! Velocity in the X, Y and Z directions:
  !To find the velocity we must first obtain the eulerian cell
  !closest to the particle.  We store these cell indices in
  !the array "cell", and the x coordinate in "xcell".  
  
  !Use faces for X, but centers for y and z.
  !When we do Y we want faces for Y and centers for X and Z.
  do axis = IAXIS, KAXIS
     if (axis .eq. IAXIS) then
        facecent = (/ FACES, CENTER, CENTER/)
        faces = (/FACEX,CENTER,CENTER/)
        particleposn(IAXIS) = particleData(POSX_PART_PROP)
     elseif (axis .eq. JAXIS) then
        facecent = (/ CENTER, FACES, CENTER/)
        faces = (/CENTER,FACEY,CENTER/)
        particleposn(JAXIS) = particleData(POSY_PART_PROP)
     elseif (axis .eq. JAXIS) then
        facecent = (/ CENTER, CENTER, FACES/)
        faces = (/CENTER,CENTER,FACEZ/)
        particleposn(KAXIS) = particleData(POSZ_PART_PROP)
     endif
     call Grid_getBlkIndexLimits(blkID, blkLimits, blkLimitsGC, faces(axis))
     coordSize = blkLimits(HIGH,axis)-blkLimits(LOW,axis)+1
     allocate(coord(coordSize), dist(coordSize))
     
     call Grid_getCellCoords(axis, blkID, facecent(axis), .false., coord, coordSize)
     
     dist = abs(coord(:) - particleposn(axis))
     cellval(axis) = minval(dist,1)
     cell(axis) = minloc(dist,1) + blkLimits(LOW,axis) - 1
           
     if (axis .eq. IAXIS) xcell = coord(minloc(dist,1))
     if (axis .eq. JAXIS) ycell = coord(minloc(dist,1))
     if (axis .eq. KAXIS) zcell = coord(minloc(dist,1))
     deallocate(coord, dist)
  end do

  !print *, "Particle:", particleposn, "in cell:", cell
  ! Now get some value for interpolated velocities in particle
!  call Grid_getBlkPtr(blkID, fcx, FACEX)
!  call Grid_getDeltas(blkID,del)
  !gr_sbBodyInfo(b) % particles(VELX_PART_PROP,localParticleCount) =  1./3.* &
  !(fcx(VELC_FACE_VAR,cell(IAXIS),cell(JAXIS),cell(KAXIS))  + & 
  !fcx(VELC_FACE_VAR,cell(IAXIS)-1,cell(JAXIS),cell(KAXIS)) + &
  !fcx(VELC_FACE_VAR,cell(IAXIS)+1,cell(JAXIS),cell(KAXIS)) )
  
  
  ! This is just one test for the interpolation from eulerian
  ! grid to Particle. ei = a variable used to find location of
  ! particle in x, in terms of what we wrote for fcx in
  ! Simulation_initBlock (x positions on VELC_FACE_VAR).
!  ei = (gr_sbBodyInfo(b) % particles(POSX_PART_PROP,localParticleCount) - xcell)/del(IAXIS) 

  ! Interpolate fcx variable with 3 points in x direction:
!  gr_sbBodyInfo(b) % particles(VELX_PART_PROP,localParticleCount) =  1./1.* &
!       (.5*ei*(ei-1)*fcx(VELC_FACE_VAR,cell(IAXIS)-1,cell(JAXIS),cell(KAXIS))  + & 
!       (1.-ei**2)*fcx(VELC_FACE_VAR,cell(IAXIS),cell(JAXIS),cell(KAXIS)) + &
!       .5*ei*(ei+1)*fcx(VELC_FACE_VAR,cell(IAXIS)+1,cell(JAXIS),cell(KAXIS)) )

  ! Similarly for y and z
!  call Grid_getBlkPtr(blkID, fcy, FACEY)
!  call Grid_getDeltas(blkID,del)
!  ej = (gr_sbBodyInfo(b) % particles(POSY_PART_PROP,localParticleCount) - ycell)/del(JAXIS)
!  gr_sbBodyInfo(b) % particles(VELY_PART_PROP,localParticleCount) =  1./1.* &
!       (.5*ej*(ej-1)*fcy(VELC_FACE_VAR,cell(IAXIS)-1,cell(JAXIS)-1,cell(KAXIS))  + &
!       (1.-ej**2)*fcy(VELC_FACE_VAR,cell(IAXIS),cell(JAXIS),cell(KAXIS)) + &
!       .5*ej*(ej+1)*fcy(VELC_FACE_VAR,cell(IAXIS),cell(JAXIS)+1,cell(KAXIS)) )

!  call Grid_getBlkPtr(blkID, fcz, FACEZ)
!  call Grid_getDeltas(blkID,del)
!  ek = (gr_sbBodyInfo(b) % particles(POSZ_PART_PROP,localParticleCount) - zcell)/del(KAXIS)
!  gr_sbBodyInfo(b) % particles(VELZ_PART_PROP,localParticleCount) =  1./1.* &
!       (.5*ek*(ek-1)*fcz(VELC_FACE_VAR,cell(IAXIS),cell(JAXIS),cell(KAXIS)-1)  + &
!       (1.-ek**2)*fcz(VELC_FACE_VAR,cell(IAXIS),cell(JAXIS),cell(KAXIS)) + &
!       .5*ek*(ek+1)*fcz(VELC_FACE_VAR,cell(IAXIS),cell(JAXIS),cell(KAXIS)+1) )

  ! Check that the interpolated variable in fcx is the
  ! actual position in x of particle.
!  write(6,'(a,i5,a,i4,a,e12.2,a,e12.2)') "Body", b, " Particle", localParticleCount, &
!       " has position", &
!       gr_sbBodyInfo(b) % particles(POSX_PART_PROP,localParticleCount), &
!       " and interpolated velocity", &
!       gr_sbBodyInfo(b) % particles(VELX_PART_PROP,localParticleCount)

!  write(6,'(a,i5,a,i4,a,e12.2,a,e12.2)') "Body", b, " Particle", localParticleCount, &
!       " has position", &
!       gr_sbBodyInfo(b) % particles(POSY_PART_PROP,localParticleCount), &
!       " and interpolated velocity", &
!       gr_sbBodyInfo(b) % particles(VELY_PART_PROP,localParticleCount)

!  write(6,'(a,i5,a,i4,a,e12.2,a,e12.2)') "Body", b, " Particle", localParticleCount, &
!       " has position", &
!       gr_sbBodyInfo(b) % particles(POSZ_PART_PROP,localParticleCount), &
!       " and interpolated velocity", &
!       gr_sbBodyInfo(b) % particles(VELZ_PART_PROP,localParticleCount)
  
  !REPRESENT LINEAR VARIATION IN X WITH 3 POINTS (AND QUADRATIC INTERPOLATION)
  
  ! Up to here interpolation of Eulerian (FACEVAR) velocities to
  ! particles velocities in VELi_PART_PROP -----------------------
  ! End of routine interpolate to particle.
  
  ! From now on, for the particle, we extrapolate a force back to Eulerian grid
  ! (FACEVAR) in IBFC_FACE_VAR
  
  ! ----- will go in routine extrapolate to eulerian stencil ---
  
!  fcx(IBFC_FACE_VAR,cell(IAXIS)-1,cell(JAXIS),cell(KAXIS)) = &
!       fcx(IBFC_FACE_VAR,cell(IAXIS)-1,cell(JAXIS),cell(KAXIS)) + &
!       gr_sbBodyInfo(b) % particles(VELX_PART_PROP,localParticleCount)
        
!  fcx(IBFC_FACE_VAR,cell(IAXIS),cell(JAXIS),cell(KAXIS)) = &
!       fcx(IBFC_FACE_VAR,cell(IAXIS),cell(JAXIS),cell(KAXIS)) + &
!       gr_sbBodyInfo(b) % particles(VELX_PART_PROP,localParticleCount) 

!  fcx(IBFC_FACE_VAR,cell(IAXIS)+1,cell(JAXIS),cell(KAXIS)) = &
!       fcx(IBFC_FACE_VAR,cell(IAXIS)+1,cell(JAXIS),cell(KAXIS)) + &
!       gr_sbBodyInfo(b) % particles(VELX_PART_PROP,localParticleCount)
        
  ! --- end of routine extrapolate to eulerian.
!  call Grid_releaseBlkPtr(blkID,fcx)
!  call Grid_releaseBlkPtr(blkID,fcy)
!  call Grid_releaseBlkPtr(blkID,fcz)
  
  return
end subroutine Grid_updateSolidBodyForces
