!!****if* source/Simulation/SimulationMain/unitTest/particleTests/VParticles/pt_initPositions
!!
!! NAME
!!    pt_initPositions
!!
!! SYNOPSIS
!!
!!    call pt_initPositions(integer(in)  :: blockID,
!!                          logical(out) :: success)
!!
!! DESCRIPTION
!!
!!    Initializes particle locations for one block in the grid.
!!
!! ARGUMENTS
!!
!!  blockID:        local block ID containing particles to create
!!
!!  success:        returns .TRUE. if positions for all particles
!!                  that should be assigned to this block have been
!!                  successfully initialized.
!!
!!***


subroutine pt_initPositions (blockID,success)
  
  use Grid_interface, ONLY : Grid_getBlkBoundBox
  use Simulation_data,ONLY : sim_initPos
  use Particles_data, ONLY:  pt_numLocal, particles, pt_maxPerProc, &
       pt_xmin, pt_ymin, pt_zmin,pt_xmax, pt_ymax, pt_zmax,&
       pt_posAttrib,pt_velNumAttrib, pt_velAttrib,pt_typeInfo, pt_meshMe
  

  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "Particles.h"
  
  integer, INTENT(in) :: blockID
  logical,intent(OUT) :: success

! ------------- local variables ------------
  real                :: xpos, ypos, zpos, rpos, bxl, byl, bzl, bxu, byu, bzu
  integer :: parts_given,p
  real, dimension(2,MDIM):: boundBox
  integer :: part_props=NPART_PROPS,k,i,j,mapType
  logical :: IsInBlock
  !-----------------------------------------
  
  p = pt_numLocal
  parts_given=CONSTANT_ONE
   

  ! Get grid geometry for this block   
  call Grid_getBlkBoundBox(blockID,boundBox)
  bxl = boundBox(LOW,IAXIS)
  bxu = boundBox(HIGH,IAXIS)
  
  if (NDIM >= 2) then
     byl = boundBox(LOW,JAXIS)
     byu = boundBox(HIGH,JAXIS)
  endif
  if (NDIM == 3) then
     bzl = boundBox(LOW,KAXIS)
     bzu = boundBox(HIGH,KAXIS)
  endif
  
  xpos=0.0;
  ypos=0.0
  zpos=0.0
  loop_x:  do i = 1, parts_given

     !xpos = (i-0.5)*dxParticle(IAXIS) + pt_initialXMin
     xpos = sim_initPos(IAXIS);
     IsInBlock = (xpos >= bxl) .and. (xpos < bxu)
          
     if (.not. IsInBlock) cycle loop_x   !! skip rest of statements if not in block
     
     loop_y:  do j = 1, parts_given
        if (NDIM >= 2) then
           ypos = sim_initPos(JAXIS);
           IsInBlock = (ypos >= byl) .and. (ypos < byu)
           if (.not. IsInBlock) cycle loop_y
        endif
        
        loop_z:  do k = 1, parts_given
           if (NDIM == 3) then
              zpos = sim_initPos(KAXIS);
              IsInBlock = (zpos >= bzl) .and. (zpos < bzu)
           endif
           
           if (IsInBlock) then
              p = p + 1
              !! Check space allocation
              if (p > pt_maxPerProc) then
                 print *,' '
                 print *,'PARAMETER pt_maxPerProc is set to ',pt_maxPerProc
                 print *,'  To avoid this crash, redimension bigger in your flash.par'
                 call Driver_abortFlash &
                      ("pt_initPositionsLattice:  Exceeded max # of particles/processor!")
              endif
              !! particle is defined, set up data structure
              
              
              particles(BLK_PART_PROP,p) = real(blockID)
              particles(PROC_PART_PROP,p) = real(pt_meshMe)
#ifdef MASS_PART_PROP
              particles(MASS_PART_PROP,p) = 1.
#endif
              particles(POSX_PART_PROP,p) = xpos
              particles(POSY_PART_PROP,p)  = ypos
              particles(POSZ_PART_PROP,p)  = zpos
              
              
           endif   !! end of IsInBlock .and. IsInSphere is true

        enddo loop_z
     enddo loop_y
  enddo loop_x
  
  pt_numLocal = p
  mapType=pt_typeInfo(PART_MAPMETHOD,1)
  !call Grid_mapMeshToParticles(particles,&
  !     part_props,BLK_PART_PROP, pt_numLocal,&
  !     pt_posAttrib,pt_velNumAttrib,pt_velAttrib,mapType)
  
  
  !-----------
  success = .TRUE.
  return
  
  !----------------------------------------------------------------------
  
end subroutine pt_initPositions


