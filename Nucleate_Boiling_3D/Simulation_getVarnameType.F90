!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! File created at setup time.  DO NOT EDIT !!!
!!

subroutine Simulation_getVarnameType(varname,vartype)
implicit none 

#include "constants.h"
#include "Flash.h"

   integer, intent(out) :: vartype
   integer, intent(in) :: varname
   integer :: ctr

   !! varnames is an array of possible input strings (with possibly extra spaces)
   !! vartypes is an array of corresponding answers

   integer, dimension(1:42), parameter :: varnames = (/&
   & PTES_VAR, RTES_VAR, AAJUNK_VAR, ALPH_VAR, CURV_VAR, &
     DELP_VAR, DFUN_VAR, DUST_VAR, MDOT_VAR, MFLG_VAR, &
     MGW1_VAR, MGW2_VAR, MGW3_VAR, MGW4_VAR, MGW5_VAR, &
     MGW6_VAR, MGW7_VAR, NRMX_VAR, NRMY_VAR, NRMZ_VAR, &
     OMGM_VAR, OMGX_VAR, OMGY_VAR, OMGZ_VAR, PFUN_VAR, &
     PREO_VAR, PRES_VAR, RHST_VAR, ROLD_VAR, SIGP_VAR, &
     SMHV_VAR, SMRH_VAR, TEMP_VAR, TMIC_VAR, TNLQ_VAR, &
     TNVP_VAR, TOLD_VAR, TVIS_VAR, VELX_VAR, VELY_VAR, &
     VELZ_VAR, VISC_VAR &
   &/)

   integer, dimension(1:42), parameter :: vartypes = (/&
   & VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC,&
     VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC,&
     VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC,&
     VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC,&
     VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC,&
     VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC,&
     VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC,&
     VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC,&
     VARTYPE_GENERIC, VARTYPE_GENERIC &
   &/)

   ! Species and mass scalars are kept in primitive form. - KW
   if (varname .ge. SPECIES_BEGIN .and. varname .le. MASS_SCALARS_END) then
      vartype = VARTYPE_PER_MASS
      return
   end if

   vartype = VARTYPE_ERROR  ! By default vartype is an ERROR
   do ctr = LBOUND(varnames,1), UBOUND(varnames,1)
      if ( varnames(ctr) .eq. varname ) vartype = vartypes(ctr)
   end do

end subroutine Simulation_getVarnameType

