!!****if* source/Simulation/SimulationMain/Nuc2Grid/sim_outputGridData
!!
!! NAME
!!  sim_outputGridData
!!
!! SYNOPSIS
!!
!!  use sim_outputGridData
!!
!! DESCRIPTION
!!
!!  Store the data for the output grid of the
!!  nucOutput-to-radTrans converter.
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!
!!   
!!
!!***

module sim_outputGridData
#include "Flash.h"
#include "constants.h"
  
  implicit none

  !! *** Runtime Parameters *** !!

  character(len=80),save :: sim_ogRadTranDataFileName
  

  integer,save :: ndimOg !dimensionality of output grid
  integer,save :: nIOg, nJOg, nKOg !grid size of output grid
  integer,save :: iOgB,iOgE, jOgB,jOgE, kOgB,kOgE !index ranges for output grid
  integer,save :: kOg2d, kOg3d
  integer,save :: iOgB0,iOgE0, jOgB0,jOgE0, kOgB0,kOgE0 !index ranges, including rim cells

  integer, parameter :: nvarsOg = NUNK_VARS !max number of vars that output grid can hold
  integer, parameter :: nvarsOgOut = 46

  !NUNK_VARS - 14 !number of output grid vars to actually write

  !When nvarsOgOut = NUNK_VARS - 13
  !velz,gaus,gpot,RPV{1,2,3},flam,gam{c,e},eint,grac,_numc_,_nup{0,1}_,pden,entr

  !When nvarsOgOut = NUNK_VARS - 14
  !velz,gaus,gpot,RPV{1,2,3},flam,gam{c,e},eint,grac,_numc_,_nup{0,1}_,pden,entr,star

  !When nvarsOgOut = NUNK_VARS - 15
  !velz,gaus,gpot,RPV{1,2,3},flam,gam{c,e},eint,grac,_numc_,_nup{0,1}_,pden,entr,file,fils

  real,allocatable,save :: outGridData(:,:,:,:) !size is (nIOg, nJOg, nKOg, nvarsOg)
  real,allocatable,save :: outLineData(:) !size is nvarsOg
  integer,allocatable,save :: outGridNumCount(:,:,:) !size is (nIOg, nJOg, nKOg)
  real,allocatable,save :: outGridMappedVol(:,:,:) !size is (nIOg, nJOg, nKOg)

  real,save :: xOgMin,xOgMax,xOgStep,xOgFact
  real,save :: yOgMin,yOgMax,yOgStep,yOgFact
  real,save :: zOgMin,zOgMax,zOgStep,zOgFact

  character(len=MAX_STRING_LENGTH),save :: geometryOgStr
  integer,save :: geometryOg !geometry of output grid

end module sim_outputGridData
