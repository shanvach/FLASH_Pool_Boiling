!!****if* source/Simulation/SimulationMain/INavierStokes/RuntimeParameters_mapStrToInt
!!
!! NAME
!!
!!  RuntimeParameters_mapStrToInt
!!
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_mapStrToInt(character(in) :: inputString(:),
!!                                     integer(out)  :: constKey)
!!
!!
!! DESCRIPTION
!!
!!  Convert a string parameter into the corresponding integer constant.
!!  The strings are defined in Config files and  provided by the flash.par file.
!!  The integer constants are defined in the header file constants.h
!!
!!  This routine is often used when mapping boundary conditions or geometry
!!  type from a string given in the flash.par to a constant key which
!!  is used by the rest of the code.
!!
!! 
!! ARGUMENTS
!!   
!!  inputString - input character string 
!!  constKey -    output integer key corresponding to inputString
!!
!! EXAMPLE
!!
!!  !  Determine the geometry requested by the flash.par
!!  call RuntimeParameters_get("geometry",pt_str_geometry)
!!  call RuntimeParameters_mapStrToInt(pt_str_geometry, pt_geometry)
!!
!!  if (pt_geometry == CARTESIAN) then
!!     .... code for rectangular domain
!!  else
!!     .... code for non-rectangular
!!  endif
!!
!!***

subroutine RuntimeParameters_mapStrToInt (inputString, constKey)

implicit none
#include "Flash.h"
#include "constants.h"

  character(len=*), intent(in) :: inputString
  integer, intent(inout) :: constKey

  constKey = NONEXISTENT

  select case (inputString)
     
  case ("periodic", "PERIODIC")
#ifdef PERIODIC
     constKey = PERIODIC
#endif

  case ("reflect", "reflecting", "REFLECT", "REFLECTING")
#ifdef REFLECTING
     constKey = REFLECTING
#endif

case ("axisymmetric", "axisymmetry", "AXISYMMETRIC", "AXISYMMETRY")
#ifdef AXISYMMETRIC
constKey = AXISYMMETRIC
#endif

case ("eqtsymmetric", "eqtsymmetry", "EQTSYMMETRIC", "EQTSYMMETRY")
#ifdef EQTSYMMETRIC
constKey = EQTSYMMETRIC
#endif
     
  case ("OUTFLOW", "neumann", "zero-gradient", "outflow")
#ifdef OUTFLOW
     constKey = OUTFLOW
#endif

   case ("diode", "DIODE")
#ifdef DIODE
     constKey = DIODE
#endif

  case ("DIRICHLET", "Dirichlet", "dirichlet")
#ifdef DIRICHLET
     constKey = DIRICHLET
#endif

  case("NEUMANN_INS","Neumann_ins","neumann_ins")
#ifdef NEUMANN_INS
     constKey = NEUMANN_INS
#endif

  case("OUTFLOW_INS","Outflow_ins","outflow_ins")
#ifdef OUTFLOW_INS
     constKey = OUTFLOW_INS
#endif

  case("NOSLIP_INS","Noslip_ins","noslip_ins")
#ifdef NOSLIP_INS
     constKey = NOSLIP_INS
#endif

  case("SLIP_INS","Slip_ins","slip_ins")
#ifdef SLIP_INS
     constKey = SLIP_INS
#endif

  case("INFLOW_INS","Inflow_ins","inflow_ins")
#ifdef INFLOW_INS
     constKey = INFLOW_INS
#endif

  case("MOVLID_INS","Movlid_ins","movlid_ins")
#ifdef MOVLID_INS
     constKey = MOVLID_INS
#endif

   case ("hydrostatic", "HYDROSTATIC")
#ifdef HYDROSTATIC
     constKey = HYDROSTATIC
#endif

   case ("hydrostatic+nvdiode")
#ifdef HYDROSTATIC_NVDIODE
     constKey = HYDROSTATIC_NVDIODE
#endif

   case ("hydrostatic+nvrefl")
#ifdef HYDROSTATIC_NVREFL
     constKey = HYDROSTATIC_NVREFL
#endif

   case ("hydrostatic+nvout")
#ifdef HYDROSTATIC_NVOUT
     constKey = HYDROSTATIC_NVOUT
#endif

   case ("hydrostatic+nvzero")
#ifdef HYDROSTATIC_NVOUT
     constKey = HYDROSTATIC_NVZERO
#endif

   case ("hydrostatic-f2", "hydrostatic-F2", "HYDROSTATIC-F2")
#ifdef HYDROSTATIC_F2
     constKey = HYDROSTATIC_F2
#endif

   case ("hydrostatic-f2+nvdiode", "hydrostatic-F2+nvdiode")
#ifdef HYDROSTATIC_F2_NVDIODE
     constKey = HYDROSTATIC_F2_NVDIODE
#endif

   case ("hydrostatic-f2+nvrefl", "hydrostatic-F2+nvrefl")
#ifdef HYDROSTATIC_F2_NVREFL
     constKey = HYDROSTATIC_F2_NVREFL
#endif

   case ("hydrostatic-f2+nvout", "hydrostatic-F2+nvout")
#ifdef HYDROSTATIC_F2_NVOUT
     constKey = HYDROSTATIC_F2_NVOUT
#endif

  case ("user", "user-defined","USER","USER-DEFINED")
#ifdef USER_DEFINED
     constKey = USER_DEFINED
#endif

  case ("CARTESIAN", "cartesian", "Cartesian")
#ifdef CARTESIAN
     constKey = CARTESIAN
#endif

  case ("polar", "POLAR")
#ifdef POLAR
     constKey = POLAR
#endif

  case ("cylindrical", "CYLINDRICAL", "Cylindrical")
#ifdef CYLINDRICAL
     constKey = CYLINDRICAL
#endif

  case ("spherical", "SPHERICAL", "Spherical")
#ifdef SPHERICAL
     constKey = SPHERICAL
#endif


  case("dens_ie","DENS_IE")
#ifdef MODE_DENS_EI
     constKey = MODE_DENS_EI
#endif

  case("dens_pres","DENS_PRES")
#ifdef MODE_DENS_PRES
     constKey = MODE_DENS_PRES
#endif

  case("dens_temp","DENS_TEMP")
#ifdef MODE_DENS_TEMP
     constKey = MODE_DENS_TEMP
#endif

#ifdef MODE_EOS_NOP
  case("eos_nop","EOS_NOP")
     constKey = MODE_EOS_NOP
#endif

  case("rt","RT")
#ifdef MODE_RT
     constKey = MODE_RT
#endif

  case("rp","RP")
#ifdef MODE_RP
     constKey = MODE_RP
#endif

  case("re","RE")
#ifdef MODE_RE
     constKey = MODE_RE
#endif

  case ("HYPRE_AMG", "hypre_amg")
#ifdef HYPRE_AMG
     constKey = HYPRE_AMG
#endif

  case("HYPRE_PCG","hypre_pcg")
#ifdef HYPRE_PCG
     constKey = HYPRE_PCG
#endif

  case("HYPRE_ILU","hypre_ilu")
#ifdef HYPRE_ILU
     constKey = HYPRE_ILU
#endif

  case ("HYPRE_BICGSTAB", "hypre_bicgstab")
#ifdef HYPRE_BICGSTAB
     constKey = HYPRE_BICGSTAB
#endif

  case ("HYPRE_GMRES", "hypre_gmres")
#ifdef HYPRE_GMRES
     constKey = HYPRE_GMRES
#endif

  case ("HYPRE_SPLIT", "hypre_split")
#ifdef HYPRE_SPLIT
     constKey = HYPRE_SPLIT
#endif

  case ("HYPRE_PARASAILS", "hypre_parasails")
#ifdef HYPRE_PARASAILS
     constKey = HYPRE_PARASAILS
#endif     

case ("HYPRE_HYBRID", "hypre_hybrid")
#ifdef HYPRE_SPLIT
     constKey = HYPRE_HYBRID
#endif 

  case ("HYPRE_NONE", "hypre_none")
#ifdef HYPRE_NONE
     constKey = HYPRE_NONE
#endif


  case DEFAULT
     constKey = NONEXISTENT
  end select
  return
end subroutine RuntimeParameters_mapStrToInt
