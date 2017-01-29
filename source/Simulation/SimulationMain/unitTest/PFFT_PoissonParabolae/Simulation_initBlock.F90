!!****if* source/Simulation/SimulationMain/unitTest/PFFT_PoissonParabolae/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!!  Reference:
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!  
!!
!! 
!!
!!***

!!REORDER(4): solnData

subroutine Simulation_initBlock(blockId)

!  use Simulation_data
  use Simulation_data, ONLY :sim_xMin,sim_xMax,sim_yMin,sim_yMax,sim_zMin,sim_zMax, &
       pfb_waven_x,pfb_waven_y,pfb_waven_z, pfb_alpha_x,pfb_alpha_y
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkBoundBox, Grid_getBlkCenterCoords, Grid_getDeltas

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  integer, intent(in) :: blockID
  !!$ ---------------------------------
 
  integer :: i, j, k
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) ::  blIndSize,blIndSizeGC

  real, dimension(MDIM)  :: coord,bsize
  real ::  boundBox(2,MDIM)
  real, pointer, dimension(:,:,:,:) :: solnData

  real :: Lx, Ly, Lz, xi, yi, zi, Phi_ijk, F_ijk

  real :: del(3)

  real :: sineI, sineJ, signI, signJ
  real :: dispI, dispJ, solfI, solfJ

  !----------------------------------------------------------------------
  Lx = sim_xMax - sim_xMin
  Ly = sim_yMax - sim_yMin  
  Lz = sim_zMax - sim_zMin


  ! Get blocks dx, dy ,dz:
  call Grid_getDeltas(blockID,del)

  ! Get Coord and Bsize for the block:
  ! Bounding box:
  call Grid_getBlkBoundBox(blockId,boundBox)
  bsize(:) = boundBox(2,:) - boundBox(1,:)

  call Grid_getBlkCenterCoords(blockId,coord)

  ! Point to Blocks centered variables:
  call Grid_getBlkPtr(blockID,solnData,CENTER)

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)


  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

           xi = coord(IAXIS) - 0.5*bsize(IAXIS) + &
                real(i - NGUARD - 1)*del(IAXIS) + 0.5*del(IAXIS)

           yi = coord(JAXIS) - 0.5*bsize(JAXIS) + &
                real(j - NGUARD - 1)*del(JAXIS) + 0.5*del(JAXIS)

           zi = coord(KAXIS) - 0.5*bsize(KAXIS) + &
                real(k - NGUARD - 1)*del(KAXIS) + 0.5*del(KAXIS)


           sineI = sin(2.*PI*xi*pfb_waven_x/Lx)
           signI = sign(1.0,sineI)
           sineJ = sin(2.*PI*yi*pfb_waven_y/Ly)
           signJ = sign(1.0,sineJ)
!!$           F_ijk  = -4.*PI**2 * ( (pfb_waven_x/Lx)**2. + (pfb_waven_y/Ly)**2. + (pfb_waven_z/Lz)**2. ) * Phi_ijk
           F_ijk  = - (signI + signJ)
           
!!$           Phi_ijk = cos(2.*PI*xi*pfb_waven_x/Lx + pfb_alpha_x) * &
!!$                     sin(2.*PI*yi*pfb_waven_y/Ly)*cos(2.*PI*zi*pfb_waven_z/Lz)

           dispI = 2*pfb_waven_x*xi/Lx - 0.5 - int(2*xi*pfb_waven_x/Lx)
           solfI = (1+signI)*(0.25-dispI*dispI) + (1-signI)*(dispI*dispI-0.25)
           solfI = solfI * Lx*Lx / (pfb_waven_x*pfb_waven_x)
           dispJ = 2*pfb_waven_y*yi/Ly - 0.5 - int(2*yi*pfb_waven_y/Ly)
           solfJ = (1+signJ)*(0.25-dispJ*dispJ) + (1-signJ)*(dispJ*dispJ-0.25)
           solfJ = solfJ * Ly*Ly / (pfb_waven_y*pfb_waven_y)

           Phi_ijk = 0.0625 * (solfI + solfJ)

           solnData(ASOL_VAR,i,j,k) = Phi_ijk

           solnData(DENS_VAR,i,j,k) = F_ijk

        enddo
     enddo
  enddo


  ! set values for u,v velocities and pressure
  solnData(DIFF_VAR,:,:,:) = 0.0
  solnData(PFFT_VAR,:,:,:) = 0.0


!!$  write(*,*) 'BlockID=',blockID
!!$  write(*,*) 'Center coordinates=',coord
!!$  write(*,*) 'Size =',bsize
 

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  return

end subroutine Simulation_initBlock
