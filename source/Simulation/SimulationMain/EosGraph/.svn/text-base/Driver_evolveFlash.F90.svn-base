!!****if* source/Simulation/SimulationMain/EosGraph/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  Driver_evolveFlash()
!!
!! DESCRIPTION
!!
!! This routine implements the time advancement scheme being used in
!! a simulation. The full implementation currently provided with FLASH is
!! a Strang splitting scheme. A single step in this driver
!! includes two sweeps, the first one in order XYZ, and
!! the second one in order ZYX. This driver works with directionally
!! split operators only. 
!!
!!  
!!***

!! Additional NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in FORTRAN
!! module Driver_data (in file Driver_data.F90). The other variables
!! are local to the specific routine and do not have the prefix "dr_"


#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

subroutine Driver_evolveFlash()

  use Eos_interface
#include "Eos.h"
#include "constants.h"
  implicit none
  integer :: i,j,k, outLUnitE, outLUnitP
  real,parameter :: xOgMin = 1e4
  real,parameter :: xOgMax = 1e11
  integer,parameter :: nIOg = 70
  integer,parameter :: eosMode = MODE_DENS_TEMP_EQUI
  real :: eosData(EOS_NUM)
  logical :: eosMask(EOS_VARS+1:EOS_NUM)
  real,allocatable :: xOgFaces(:)
  real :: xOgStep, xOgFact
  real :: density,temp,x,y,y1,y2,y3

  eosMask(:) = .FALSE.
  eosMask(EOS_EINTION) = .TRUE.
  eosMask(EOS_EINTELE) = .TRUE.
  eosMask(EOS_EINTRAD) = .TRUE.
  eosMask(EOS_PRESION) = .TRUE.
  eosMask(EOS_PRESELE) = .TRUE.
  eosMask(EOS_PRESRAD) = .TRUE.

  allocate(xOgFaces(0:nIOg+2))

  if (xOgMin .LE. 0.0) then     !linear spacing...

     xOgStep = (xOgMax - xOgMin) / nIOg
     do i=0,nIOg+2
        xOgFaces(i) = xOgMin + (real(i)-1.0) * xOgStep 
     end do

  else                          !log spacing...
     xOgFact = 10.0**((alog10(xOgMax) - alog10(xOgMin)) / nIOg)
     do i=0,nIOg+2
        xOgFaces(i) = xOgMin * xOgFact**(real(i)-1.0)
     end do
  end if
  
  outLUnitE = 30
  open(outLUnitE,name='EOSdumpE.dat',status='UNKNOWN')
  outLUnitP = 31
  open(outLUnitP,name='EOSdumpP.dat',status='UNKNOWN')

  density = 1.0e3

  do i = 1,nIOg+1

     x = xOgFaces(i)
     eosData(EOS_DENS) = density
     eosData(EOS_TEMP) = x
     eosData(EOS_ABAR) = 1.0
     eosData(EOS_ZBAR) = 1.0
     call Eos(eosMode,1,eosData,mask=eosMask)

     y = eosData(EOS_EINT)
     y1 = eosData(EOS_EINTION)
     y2 = eosData(EOS_EINTELE)
     y3 = eosData(EOS_EINTRAD)
     write(outLUnitE,1) x,y,y1,y2,y3

     y = eosData(EOS_PRES)
     y1 = eosData(EOS_PRESION)
     y2 = eosData(EOS_PRESELE)
     y3 = eosData(EOS_PRESRAD)
     write(outLUnitP,1) x,y,y1,y2,y3
1    format(5(1x,1PG12.3))

  end do

  return
  
end subroutine Driver_evolveFlash



