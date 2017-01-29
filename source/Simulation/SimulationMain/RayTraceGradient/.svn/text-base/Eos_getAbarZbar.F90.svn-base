!!****if* source/Simulation/SimulationMain/RayTraceGradient/Eos_getAbarZbar
!! NAME
!!
!!  Eos_getAbarZbar
!! 
!! SYNOPSIS
!!
!!  call Eos_getAbarZbar(real(IN)  :: solnVec(NUNK_VARS),
!!                     integer(IN) :: pos(MDIM),
!!                     integer(IN) :: vecLen,
!!                  real, pointer  :: solnData(:,:,:,:),
!!                     integer(IN) :: gridDataStruct,
!!                     real(OUT)   :: eosData(:),
!!             optional,real(OUT)  :: abar,
!!             optional,real(OUT)  :: zbar,
!!             optional,real(OUT)  :: sumY,
!!             optional,real(OUT)  :: Ye,
!!             optional,real(IN)   :: massFrac(:) )
!!
!!
!!
!! DESCRIPTION
!!
!! Eos_getAbarZbar gets Abar, Zbar, SumY, and/or Ye for one cell.
!!
!!
!!  ARGUMENTS 
!!
!!   
!!   solnVec : the solution vector for one cell
!!   massFrac : this is an optional argument which may be used when there is more 
!!              than one species in the simulation
!!
!!  EXAMPLE 
!!
!!  NOTES
!!
!!
!!      This interface is defined in Fortran Module 
!!      Eos_interface. All functions calling this routine should include
!!      a statement like
!!      use Eos_interface, ONLY : Eos_putData
!!
!!
!!  SEE ALSO
!!
!!     Eos
!!
!!
!!***

#include "Flash.h"


subroutine Eos_getAbarZbar(solnVec,abar,zbar,sumY,Ye,massFrac)

  implicit none
  
!#include "Eos.h"
!#include "Eos_map.h"
!#include "constants.h"
  
  real, OPTIONAL,dimension(NUNK_VARS),intent(IN) :: solnVec
  real, OPTIONAL,                    intent(OUT) :: abar, zbar, Ye, sumY
  real, OPTIONAL,dimension(NSPECIES), intent(IN) :: massFrac


  if (present(abar)) abar = 1.0
  if (present(zbar)) zbar = 1.0
  if (present(sumY)) sumY = 0.0
  if (present(Ye))   Ye   = 0.0

  return
end subroutine Eos_getAbarZbar



