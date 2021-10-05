!!****if* source/physics/SolidMechanics/SolidMechanicsMain/SurfaceInteraction/sm_surf_MapParticles
!!
!! NAME
!!
!!
!!
!! SYNOPSIS
!!
!!  
!! VARIABLES
!!
!!
!! DESCRIPTION
!! 
!!
!!***
#include "constants.h"
#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_surf_MapParticles(ibd)
  use SolidMechanics_data, only: sm_bodyInfo
  use sm_element_interface, only: sm_el02_mapParticles, sm_el10_mapParticles, &
                                  sm_el01_mapParticles
  use gr_sbData, ONLY : gr_sbBodyInfo
  use Driver_interface, ONLY : Driver_abortFlash
  use Simulation_data, only: sim_ibm_x, sim_ibm_y, sim_ibm_z
  use Driver_data, only: dr_simTime
  implicit none
  
  ! IO variables
  integer, intent(in) :: ibd

  ! Internal Variables
  integer :: numPart, e, ptelem, max_ptelem, nel, p
  real, allocatable, dimension(:) :: xpos,ypos,zpos,xvel,yvel,zvel,areai
  real, allocatable, dimension(:) :: xacc,yacc,zacc,xnrm,ynrm,znrm
  integer, allocatable, dimension(:) :: loc_num
  real :: x1,y1,x2,y2,x3,y3
  real :: z1,z2,z3

  ! Maximum number of particles per element:
  nel = sm_bodyInfo(ibd)%ws_nel
  max_ptelem = maxval(sm_bodyInfo(ibd)%ws_ptelem(1:nel))

  allocate( xpos(max_ptelem), xvel(max_ptelem), xacc(max_ptelem), xnrm(max_ptelem) )
  allocate( ypos(max_ptelem), yvel(max_ptelem), yacc(max_ptelem), ynrm(max_ptelem) )
  allocate( areai(max_ptelem) )
  allocate( loc_num(max_ptelem) )
#if NDIM == MDIM
  allocate( zpos(max_ptelem), zvel(max_ptelem), zacc(max_ptelem), znrm(max_ptelem) )
#endif
  numPart = 0
  
  do e=1,nel

     ptelem = sm_BodyInfo(ibd)%ws_ptelem(e)

     select case( sm_bodyInfo(ibd)%ws_eltype(e) )

#if NDIM == MDIM

        case( THREE_NODE_TRIANGLE )
           call sm_el02_mapParticles(sm_bodyInfo(ibd), e, ptelem, &
                                     xpos(1:ptelem),ypos(1:ptelem),zpos(1:ptelem), &
                                     xvel(1:ptelem),yvel(1:ptelem),zvel(1:ptelem), &
                                     xacc(1:ptelem),yacc(1:ptelem),zacc(1:ptelem), &
                                     xnrm(1:ptelem),ynrm(1:ptelem),znrm(1:ptelem), &
                                     areai(1:ptelem), loc_num(1:ptelem) )

        case( NINE_NODE_QUADRILATERAL )
           call sm_el10_mapParticles(sm_bodyInfo(ibd), e, ptelem, &
                                     xpos(1:ptelem),ypos(1:ptelem),zpos(1:ptelem), &
                                     xvel(1:ptelem),yvel(1:ptelem),zvel(1:ptelem), &
                                     xacc(1:ptelem),yacc(1:ptelem),zacc(1:ptelem), &
                                     xnrm(1:ptelem),ynrm(1:ptelem),znrm(1:ptelem), &
                                     areai(1:ptelem), loc_num(1:ptelem) )

#else /* 2D */

        case( TWO_NODE_LINE )
           call sm_el01_mapParticles(sm_bodyInfo(ibd), e, ptelem, &
                                   xpos(1:ptelem),ypos(1:ptelem), &
                                   xvel(1:ptelem),yvel(1:ptelem), &
                                   xacc(1:ptelem),yacc(1:ptelem), &
                                   xnrm(1:ptelem),ynrm(1:ptelem), &
                                   areai(1:ptelem), loc_num(1:ptelem) )


#endif

        case default
           call Driver_abortFlash("Element type not recognized.")

     end select
  
     ! Load elements Particles data into particles array:
     gr_sbBodyInfo(ibd) % particles(BLK_PART_PROP,numPart+1:numPart+ptelem)  = UNKNOWN
     gr_sbBodyInfo(ibd) % particles(POSX_PART_PROP,numPart+1:numPart+ptelem) = xpos(1:ptelem)
     gr_sbBodyInfo(ibd) % particles(POSY_PART_PROP,numPart+1:numPart+ptelem) = ypos(1:ptelem)
     gr_sbBodyInfo(ibd) % particles(VELX_PART_PROP,numPart+1:numPart+ptelem) = xvel(1:ptelem)
     gr_sbBodyInfo(ibd) % particles(VELY_PART_PROP,numPart+1:numPart+ptelem) = yvel(1:ptelem)
     gr_sbBodyInfo(ibd) % particles(ACCX_PART_PROP,numPart+1:numPart+ptelem) = xacc(1:ptelem)
     gr_sbBodyInfo(ibd) % particles(ACCY_PART_PROP,numPart+1:numPart+ptelem) = yacc(1:ptelem)
     gr_sbBodyInfo(ibd) % particles(NMLX_PART_PROP,numPart+1:numPart+ptelem) = xnrm(1:ptelem)
     gr_sbBodyInfo(ibd) % particles(NMLY_PART_PROP,numPart+1:numPart+ptelem) = ynrm(1:ptelem)

     gr_sbBodyInfo(ibd) % particles(PRES_PART_PROP,numPart+1:numPart+ptelem) = 0.0
     gr_sbBodyInfo(ibd) % particles(FXVI_PART_PROP,numPart+1:numPart+ptelem) = 0.0
     gr_sbBodyInfo(ibd) % particles(FYVI_PART_PROP,numPart+1:numPart+ptelem) = 0.0

#if NDIM == MDIM
     gr_sbBodyInfo(ibd) % particles(POSZ_PART_PROP,numPart+1:numPart+ptelem) = zpos(1:ptelem)
     gr_sbBodyInfo(ibd) % particles(VELZ_PART_PROP,numPart+1:numPart+ptelem) = zvel(1:ptelem)
     gr_sbBodyInfo(ibd) % particles(ACCZ_PART_PROP,numPart+1:numPart+ptelem) = zacc(1:ptelem)
     gr_sbBodyInfo(ibd) % particles(NMLZ_PART_PROP,numPart+1:numPart+ptelem) = znrm(1:ptelem)
     gr_sbBodyInfo(ibd) % particles(FZVI_PART_PROP,numPart+1:numPart+ptelem) = 0.0
#endif

     ! Bookkeeping: element number, and local linear patch number:
     gr_sbBodyInfo(ibd) % particles(ELEM_PART_PROP, numPart+1:numPart+ptelem) = real( e )
     gr_sbBodyInfo(ibd) % particles(PLOC_PART_PROP, numPart+1:numPart+ptelem) = real( loc_num(1:ptelem) )

     gr_sbBodyInfo(ibd) % particles(FUL_PART_PROP,numPart+1:numPart+ptelem)  = 0.
     gr_sbBodyInfo(ibd) % particles(FVL_PART_PROP,numPart+1:numPart+ptelem)  = 0. 
#if NDIM == MDIM
     gr_sbBodyInfo(ibd) % particles(FWL_PART_PROP,numPart+1:numPart+ptelem)  = 0.
#endif
     gr_sbBodyInfo(ibd) % particles(TAG_PART_PROP,numPart+1:numPart+ptelem)  = 1.
     gr_sbBodyInfo(ibd) % particles(AREA_PART_PROP,numPart+1:numPart+ptelem) = areai(1:ptelem)

     do p = numPart+1,numPart+ptelem
        gr_sbBodyInfo(ibd) % particles(GLOB_PART_PROP,p) = p !local particle counter in solid body.
     enddo

     numPart = numPart + ptelem
     
     if(e .eq. 1) then
      x1 = xpos(1)
      y1 = ypos(1)
     end if

     if(e .eq. 15) then
      x2 = xpos(1)
      y2 = ypos(1)
     end if

     if(e .eq. 30) then
      x3 = xpos(1)
      y3 = ypos(1)
     end if

  enddo

  !sim_ibm_x = ((x2**2+y2**2)*(y3-y1) + (x1**2+y1**2)*(y2-y3) + (x3**2+y3**2)*(y1-y2)) &
  !          / (2*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)))
 
  !sim_ibm_y = ((x2**2+y2**2)*(x3-x1) + (x1**2+y1**2)*(x2-x3) + (x3**2+y3**2)*(x1-x2)) &
  !          / (2*((y2-y1)*(x3-x1)-(y3-y1)*(x2-x1)))

  !sim_ibm_x = xpos(1) - 0.5
  !sim_ibm_y = ypos(1)
 
  sim_ibm_x = 0.0
  sim_ibm_y = -1.0*dr_simTime
  sim_ibm_z = 0.0

  print *,"P1: ", x1, y1, z1
  print *,"P2: ", x2, y2, z2
  print *,"P3: ", x3, y3, z3

  print *,"Center: ",sim_ibm_x, sim_ibm_y, sim_ibm_z

  ! deallocate temp storage
  deallocate( xpos, xvel, xacc, xnrm )
  deallocate( ypos, yvel, yacc, ynrm )
  deallocate( areai , loc_num )
#if NDIM == MDIM
  deallocate( zpos, zvel, zacc, znrm )
#endif
  
  return

end subroutine sm_surf_MapParticles
