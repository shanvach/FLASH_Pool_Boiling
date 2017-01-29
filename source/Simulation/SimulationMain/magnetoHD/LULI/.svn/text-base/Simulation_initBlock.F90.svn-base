!!****if* source/Simulation/SimulationMain/magnetoHD/LULI/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!! 
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!  
!!
!!
!!***

subroutine Simulation_initBlock(blockId)
  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_putPointData,      &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getDeltas
  use Driver_interface, ONLY: Driver_abortFlash
  use RadTrans_interface, ONLY: RadTrans_mgdEFromT

  implicit none

#include "constants.h"
#include "Flash.h"

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)

  integer, intent(in) :: blockId
  
  integer :: i, j, k
  integer :: blkLimits(2, MDIM)
  integer :: blkLimitsGC(2, MDIM)
  real, allocatable :: xcent(:), ycent(:), zcent(:)
  real :: tradActual
  real :: rho, tele, trad, tion
  integer :: species
  real, pointer, dimension(:,:,:,:) :: solnData, facexData, faceyData

  integer :: ndx, ndy, ndz, nsctot
  integer :: nsct, nscc, l, m, n
  real :: xsc, ysc, zsc, targFrac, chamFrac, norm
  real :: delta(MDIM)
  real :: r, phi

  ndx = sim_ndiv
  ndy = 1
  ndz = 1
  if(NDIM >= 2) ndy = sim_ndiv
  if(NDIM == 3) ndz = sim_ndiv
  nsctot = ndx * ndy * ndz

  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  allocate(xcent(blkLimitsGC(HIGH, IAXIS)))
  call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., &
       xcent, blkLimitsGC(HIGH, IAXIS))
  allocate(ycent(blkLimitsGC(HIGH, JAXIS)))
  call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., &
       ycent, blkLimitsGC(HIGH, JAXIS))
  allocate(zcent(blkLimitsGC(HIGH, KAXIS)))
  call Grid_getCellCoords(KAXIS, blockId, CENTER, .true., &
       zcent, blkLimitsGC(HIGH, KAXIS))

  call Grid_getDeltas(blockId, delta)

  !------------------------------------------------------------------------------

  call Grid_getBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
  endif
#endif

  ! Loop over cells and set the initial state
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           nsct = 0
           do l = 1, ndx
              do m = 1, ndy
                 do n = 1, ndz
                    xsc = (xcent(i)-delta(IAXIS)/2) + (delta(IAXIS)*(l-1))/ndx
                    ysc = (ycent(j)-delta(JAXIS)/2) + (delta(JAXIS)*(m-1))/ndy
                    zsc = (zcent(k)-delta(KAXIS)/2) + (delta(KAXIS)*(n-1))/ndz
                    call findSpec(xsc, ysc, zcent(k), species)
                    if(species == TARG_SPEC) nsct = nsct + 1
                 end do
              end do
           end do

           nscc = nsctot - nsct

           ! Average Cell Density:                                                               
           rho = (sim_rhoTarg * nsct + sim_rhoCham * nscc)/nsctot

           ! Mass Fractions:                                                                     
           targFrac = (sim_rhoTarg*nsct)/(rho * nsctot)
           chamFrac = 1.0 - targFrac

           ! Set temperatures:
           tele = (sim_teleTarg * nsct + sim_teleCham * nscc)/nsctot
           tion = (sim_tionTarg * nsct + sim_tionCham * nscc)/nsctot
           trad = (sim_tradTarg * nsct + sim_tradCham * nscc)/nsctot

           ! Set all initial conditions:
           solnData(DENS_VAR,i,j,k) = rho
           solnData(TEMP_VAR,i,j,k) = tele
           solnData(TION_VAR,i,j,k) = tion
           solnData(TELE_VAR,i,j,k) = tele

           targFrac = targFrac + sim_smallX
           chamFrac = chamFrac + sim_smallX
           solnData(TARG_SPEC,i,j,k) = targFrac/(targFrac + chamFrac)
           solnData(CHAM_SPEC,i,j,k) = chamFrac/(targFrac + chamFrac)

           solnData(HEAT_MSCALAR,i,j,k) = solnData(TARG_SPEC,i,j,k)
           ! if(xcent(i) < 0.0125) then
           !    solnData(HEAT_MSCALAR,i,j,k) = solnData(TARG_SPEC,i,j,k)
           ! else
           !    solnData(HEAT_MSCALAR,i,j,k) = 0.0
           ! end if

#ifdef BDRY_VAR
           solnData(BDRY_VAR,i,j,k) = -1.0
           if((sim_useMesh .eqv. .true.) .and. (trim(sim_meshGeom) == "flat")) then
              if(ycent(j) >= 1.0 .and. ycent(j) <= 1.02) then
                 if(  (xcent(i) > 0.0 .and. xcent(i) < 0.1) .or. &
                      (xcent(i) > 0.2 .and. xcent(i) < 0.3) .or. &
                      (xcent(i) > 0.4 .and. xcent(i) < 0.5) .or. &
                      (xcent(i) > 0.6 .and. xcent(i) < 0.7) .or. &
                      (xcent(i) > 0.8 .and. xcent(i) < 0.9) .or. & 
                      (xcent(i) > 1.0 .and. xcent(i) < 1.1) .or. & 
                      (xcent(i) > 1.2 .and. xcent(i) < 1.3) .or. & 
                      (xcent(i) > 1.2 .and. xcent(i) < 1.3) .or. & 
                      (xcent(i) > 1.4 .and. xcent(i) < 1.5) .or. & 
                      (xcent(i) > 1.6 .and. xcent(i) < 1.7) .or. & 
                      (xcent(i) > 1.8 .and. xcent(i) < 1.9) ) then
                    solnData(BDRY_VAR,i,j,k) = 1.0
                 end if
              end if
           end if

           if((sim_useMesh .eqv. .true.) .and. (trim(sim_meshGeom) == "curved")) then
              r = sqrt(xcent(i)**2 + ycent(j)**2)
              phi = atan(ycent(j)/xcent(i)) * 180.0/PI

              if(r >= 1.0 .and. r <= 1.02) then
                 if(  (phi >  0 .and. phi <  6) .or. &
                      (phi > 12 .and. phi < 18) .or. &
                      (phi > 24 .and. phi < 30) .or. &
                      (phi > 36 .and. phi < 42) .or. &
                      (phi > 48 .and. phi < 54) .or. &
                      (phi > 60 .and. phi < 66) .or. &
                      (phi > 72 .and. phi < 78) .or. &
                      (phi > 84 .and. phi < 90) ) then
                    solnData(BDRY_VAR,i,j,k) = 1.0
                 end if
              end if
           end if

#endif

           ! Set up radiation energy density:
           call RadTrans_mgdEFromT(blockId, (/i,j,k/), trad, tradActual)
           solnData(TRAD_VAR,i,j,k) = tradActual

           !! magnetic fields
#ifdef MAGX_VAR
           solnData(MAGX_VAR,i,j,k) = sim_Bx
#endif
#ifdef MAGY_VAR
           solnData(MAGY_VAR,i,j,k) = sim_By
#endif
#ifdef MAGZ_VAR
           solnData(MAGZ_VAR,i,j,k) = sim_Bz
#endif


#if NFACE_VARS > 0
           if (sim_killdivb) then
              facexData(MAG_FACE_VAR,i,j,k)= sim_Bx
              faceyData(MAG_FACE_VAR,i,j,k)= sim_By
           endif
#endif

        enddo
     enddo
  enddo


#if NFACE_VARS > 0
  if (sim_killdivb) then
     
     do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
        i=blkLimitsGC(HIGH,IAXIS)+1
        do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           facexData(MAG_FACE_VAR,i,j,k)= sim_Bx
        enddo
        
        j=blkLimitsGC(HIGH,JAXIS)+1
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           faceyData(MAG_FACE_VAR,i,j,k) = sim_By
        enddo
     enddo
     
  endif
#endif


  deallocate(xcent)
  deallocate(ycent)
  deallocate(zcent)

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
  endif
#endif

  return

contains

  subroutine findSpec(x,y,z, spec)
    implicit none

    real, intent(in) :: x
    real, intent(in) :: y
    real, intent(in) :: z
    
    integer, intent(out) :: spec
    
    spec = CHAM_SPEC

    select case (NDIM)
    case (1)

       select case (sim_geometry)
       case (SPHERICAL)
          if (sim_targetGeom == "sphere") then
             if( x < sim_targetRadius ) &
                  spec = TARG_SPEC
          else
             call Driver_abortFlash("[Simulation_initBlock] Unsupported sim_targetGeom")
          end if

       case default
          call Driver_abortFlash("[Simulation_initBlock] Unsupported geometry")
       end select

    case (2)       

       select case (sim_geometry)
       case (CARTESIAN)

          if(sim_targetGeom == "cylinder") then
             
             r = sqrt(x**2 + y**2)
             phi = atan(y/x)

             if(r <= sim_targetRadius * (1. + sim_skewFactor*cos(2*phi))) then
                spec = TARG_SPEC
             end if
          else
             call Driver_abortFlash("[Simulation_initBlock] Unsupported sim_targetGeom")
          end if

       case (CYLINDRICAL)

          
          if(sim_targetGeom == "cylinder") then
             
             if(abs(x-sim_targetOffset) <= sim_targetRadius .and. &
                  abs(y) <= sim_targetHeight ) then
                spec = TARG_SPEC
             end if

          elseif (sim_targetGeom == "sphere") then
             
             if(sqrt(x**2 + (y-sim_targetZOffset)**2) <= sim_targetRadius) then
                spec = TARG_SPEC
             end if
             
          elseif (sim_targetGeom == "2sphere") then

             if(sqrt(x**2 + (y-sim_targetZOffset)**2) <= sim_targetRadius .or. & 
                  sqrt(x**2 + (y+sim_targetZOffset)**2) <= sim_targetRadius ) then
                spec = TARG_SPEC
             end if

          else
             call Driver_abortFlash("[Simulation_initBlock] Unsupported sim_targetGeom")
          end if

       case default 
          call Driver_abortFlash("[Simulation_initBlock] Unsupported geometry")
       end select

    case (3)

       select case (sim_geometry)
       case (CARTESIAN)
          if (sim_targetGeom == "sphere") then
             if(sqrt(x**2 + y**2 + z**2) <= sim_targetRadius) &
                  spec = TARG_SPEC
          elseif (sim_targetGeom == "cylinder") then
             if(sqrt(x**2 + y**2) < sim_targetRadius .and. &
                  abs(z) <= sim_targetHeight ) &
                  spec = TARG_SPEC
          else
             call Driver_abortFlash("[Simulation_initBlock] Unsupported sim_targetGeom")
          end if

       case default 
          call Driver_abortFlash("[Simulation_initBlock] Unsupported geometry")
       end select

    end select
           
  end subroutine findSpec


end subroutine Simulation_initBlock
