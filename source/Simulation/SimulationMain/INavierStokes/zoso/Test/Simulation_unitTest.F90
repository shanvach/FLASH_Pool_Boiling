subroutine Simulation_unitTest(blockCount, blockList,timeEndAdv, dt)

#include "Flash.h"
#include "ImBound.h"
#include "SolidMechanics.h"

#ifdef FLASH_GRID_PARAMESH
  use physicaldata, ONLY : force_consistency,        &
                           interp_mask_unk_res,      &
                           interp_mask_facex_res,    &
                           interp_mask_facey_res,    &
                           interp_mask_facez_res,    &
                           interp_mask_unk,      &
                           interp_mask_facex,    &
                           interp_mask_facey,    &
                           interp_mask_facez
  use workspace, ONLY :    interp_mask_work                           
#endif    
  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
                             Grid_getListOfBlocks, &
                             Grid_getDeltas,         &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_fillGuardCells,    &
                             Grid_putFluxData,       &
                             Grid_getFluxData,       &
                             Grid_conserveFluxes,    &
                             Grid_conserveField,     &
                             Grid_updateRefinement,  &
                             Grid_solvePoisson, Grid_getBlkBoundBox, Grid_getBlkCenterCoords
  
  use gr_interface, ONLY : gr_findMean

  use ins_interface, only  :  ins_vt,&
                           ins_rhs3d,&
                           ins_rhs2d,&
                       ins_predictor,&
                      ins_divergence,&
                       ins_corrector,&
                         ins_fluxfix,&
                       ins_fluxfix_p,&
                   ins_computeQinout,&
                   ins_rescaleVelout,&
                   ins_convectVelout,&
              ins_setInterpValsGcell,&
                      ins_UstarStats,&
                  ins_pressgradients

  use IncompNS_data, ONLY : ins_isgs, ins_invRe, ins_intschm, ins_prescoeff, ins_meshMe,&
                            ins_restart, ins_nstep, ins_Qin, ins_Qout, ins_predcorrflg, &
                            ins_convvel, ins_alf, ins_gam, ins_rho, ins_gama, ins_alfa, &
                            ins_rhoa, ins_outflowgridChanged, ins_tlevel, &
                            ins_vardt, rkstep, ins_intschm_type

  use Grid_Data, ONLY : gr_domainBC 
    
  use Timers_interface, ONLY : Timers_start, Timers_stop

  use ImBound_interface, ONLY : ImBound
 
  use SolidMechanics_interface, only : SolidMechanics

  use Driver_interface, only : Driver_getNStep

  implicit none

#include "constants.h"
#include "IncompNS.h"
!#ifdef FLASH_GRID_PARAMESH
!#include "Multigrid.h"
!#endif
  include "Flash_mpi.h"

  !! ---- Argument List ----------------------------------
  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList !blockCount
  real,    INTENT(IN) :: timeEndAdv,dt
  !! -----------------------------------------------------
  
  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)

  integer :: lb,blockID,i,j,k

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(6) :: bc_types

  real :: xcell, ycell, poisfact, DIFF, ERR

  real, parameter :: Pi = 3.141592

  real, dimension(2,6) :: bc_values = 0.0
  real, dimension(MDIM)  :: coord, bsize, del
  real, dimension(2,MDIM) :: boundBox

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  ! Loop over all blocks
  do lb=1,blockCount
     blockID = blockList(lb)
     
     ! Get the cell center and block limit information for current block
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)

     ! Call for coor and del information
     call Grid_getBlkCenterCoords(blockId,coord)
     call Grid_getDeltas(blockID,del)

     ! Bounding box:
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(1:NDIM) = boundBox(2,1:NDIM) - boundBox(1,1:NDIM)

     ! Loop over all coordinates in that block
     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
            
        ! Calculate the x and y coordinates from the index for i and j
        xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                real(i - NGUARD - 1)*del(IAXIS) +   &
                0.5*del(IAXIS)

        ycell = coord(JAXIS) - bsize(JAXIS)/2.0 +   &
                real(j - NGUARD - 1)*del(JAXIS)  +  &
                0.5*del(JAXIS)

        ! Analytical
        solnData(TEMP_VAR,i,j,1)  = COS(2 * Pi * xcell) * COS(2 * Pi * ycell)

        ! Numerical
        solnData(DUST_VAR,i,j,1)  = -8 * Pi**2 * COS(2 * Pi * xcell) * & 
                                          COS(2 * Pi * ycell)  
        end do

     end do
     
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
          
  enddo

  ! Apply Boundry Conditions for DUST_VAR and TEMP_VAR
  gcMask = .FALSE.
  gcMask(DUST_VAR) = .TRUE.                           
  gcMask(TEMP_VAR) = .TRUE.                            
  
  ! Fill guard cells for DUST_VAR and TEMP_VAR 
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
    maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)             

  ! PRES_VAR = f(x,y), TEMP_VAR = phi
  call Timers_start("Grid_solvePoisson")
  ! SOLUTION OF POISSON EQUATION:
  ! -------- -- ------- -------- --- --------
  poisfact = 1.0 
  bc_types(:) = GRID_PDE_BND_NEUMANN !MG_BND_NEUMANN
  call Grid_solvePoisson(PRES_VAR, DUST_VAR, bc_types, bc_values, poisfact) 
  call Timers_stop("Grid_solvePoisson")

  ! Apply Boundry Conditions for PRES_VAR
  gcMask = .FALSE.
  gcMask(PRES_VAR) = .TRUE.                           
  
  ! Fill guard cells for PRES_VAR
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
    maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)             

  ! Loop to calculate error
  do lb=1,blockCount
     blockID = blockList(lb)
     
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)
     
     ! Calculate Error
     solnData(TVIS_VAR,:,:,1) = ABS(solnData(PRES_VAR,:,:,1) - solnData(TEMP_VAR,:,:,1))
     
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  end do

END SUBROUTINE Simulation_unitTest

