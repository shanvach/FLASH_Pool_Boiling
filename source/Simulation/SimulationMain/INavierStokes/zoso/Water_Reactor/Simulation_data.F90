!!****if* source/Simulation/SimulationMain/INavierStokes/2D/IB_Cyl_parallelIBVP_HYPRE_VD_halfDiamBel/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  Simulation_data()
!!
!! DESCRIPTION
!!
!!  Stores the local data for Simulation setup: INS-iso-turb
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"
#include "Flash.h"

  !! *** Runtime Parameters *** !!
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax
  real, save    :: sim_zMin, sim_zMax
  logical, save :: sim_gCell

  integer, save :: sim_meshMe
  real, dimension(MDIM), save :: sim_initPos

  real, save :: sim_qIn1, sim_qIn2, sim_qOut1, sim_qOut2

  real, save :: sim_waveA

  real, save, dimension(500) :: sim_nuc_site_x, sim_nuc_site_y, sim_nuc_site_z, sim_nuc_radii

  integer, save :: sim_nucSiteDens

  real, save :: sim_Tbulk

  real, save :: sim_sinkB

  real, save :: sim_xmin1, sim_xmax1, sim_xmin2, sim_xmax2
  real, save :: sim_ymin1, sim_ymax1, sim_ymin2, sim_ymax2

  real, save :: sim_vel(NXB+2*NGUARD,NZB+2*NGUARD*K3D,MAXBLOCKS)
  integer, save :: sim_vel_flg(NXB+2*NGUARD,NZB+2*NGUARD*K3D,MAXBLOCKS)
  integer, save :: sim_pres_flg(NXB+2*NGUARD,NZB+2*NGUARD*K3D,MAXBLOCKS)

  real, save :: sim_vel_top(NXB+2*NGUARD,NZB+2*NGUARD*K3D,MAXBLOCKS)
  integer, save :: sim_vel_flg_top(NXB+2*NGUARD,NZB+2*NGUARD*K3D,MAXBLOCKS)
  integer, save :: sim_pres_flg_top(NXB+2*NGUARD,NZB+2*NGUARD*K3D,MAXBLOCKS)
  real, save :: sim_dfun_top(NXB+2*NGUARD,NZB+2*NGUARD*K3D,MAXBLOCKS)

  real, save :: sim_xjet(2), sim_xjet_top(2)
 
end module Simulation_data
