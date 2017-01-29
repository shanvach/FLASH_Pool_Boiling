!!****if* source/physics/sourceTerms/Cool/CoolMain/SutherlandDopita/Cool_init
!!
!! NAME
!!
!!  Cool_init
!!
!! SYNOPSIS
!!
!!  Cool_init()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   
!!
!! AUTOGENROBODOC
!!
!!
!!***


subroutine Cool_init()

  use Cool_data, ONLY :cl_fctn_file,cl_Nmax,cl_logT, cl_ne, cl_nH, &
       cl_nt, cl_logLnet, cl_logLnorm, &
       cl_logU, cl_logtau, cl_P12, cl_rho24, cl_Ci, cl_mubar, &
       cl_dlogLdlogU, cl_dnedlogU, cl_dntdlogU,cl_N,cl_Xin, cl_Abar, useCool,&
       cl_rho, cl_gamma, cl_smalle, cl_Boltzmann, cl_AMU, cl_meshMe


  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

#include "constants.h"

  implicit none

  
  real    :: dlogUinv
  integer :: i

!==============================================================================

  ! Everybody should know these
  call Driver_getMype(MESH_COMM,cl_meshMe)

  call RuntimeParameters_get ("cool_fctn_file", cl_fctn_file)
  call PhysicalConstants_get ("proton mass", cl_AMU)
  call PhysicalConstants_get ("Boltzmann", cl_Boltzmann)
  call RuntimeParameters_get ("smalle", cl_smalle)
  call RuntimeParameters_get ("useCool",useCool)  
  if (.not. useCool) then
     write(6,*)'WARNING:  You have included the Cool unit but have set '
     write(6,*)'   the runtime parameter useCool to FALSE'
     write(6,*)'   No cooling will occur but Cool_init will continue.'
  end if

  
  ! Now allocate space for the arrays.
  
  if(useCool) then

     if (cl_meshMe .EQ. MASTER_PE) &
          print *, "init_cool:  attempting to read cooling function file '", &   ! if MASTER_PE
          cl_fctn_file(1:len_trim(cl_fctn_file)), "'..."

     call readSdTable (cl_fctn_file(1:len_trim(cl_fctn_file)), &
          cl_logT, cl_ne, cl_nH, cl_nt, cl_logLnet, cl_logLnorm, &
          cl_logU, cl_logtau, cl_P12, cl_rho24, cl_Ci, cl_mubar, cl_Nmax, CL_N)
     if (CL_N < 1) then
        print *, cl_meshMe, ":  init_cool:  could not read '", cl_fctn_file, "'!"  ! shows failure not on MASTER_PE
        call Driver_abortFlash("init_cool:  unable to read cooling function file")
     endif
     
     do i = 1, cl_N-1
        dlogUinv      = 1. / (cl_logU(i+1) - cl_logU(i))
        cl_dlogLdlogU(i) = (cl_logLnorm(i+1) - cl_logLnorm(i)) * dlogUinv
        cl_dnedlogU(i)   = (cl_ne(i+1) - cl_ne(i)) * dlogUinv
        cl_dntdlogU(i)   = (cl_nt(i+1) - cl_nt(i)) * dlogUinv
     enddo
     
     if (cl_meshMe.EQ.MASTER_PE) print *, "init_cool:  successfully read file."
     
  end if

  return
end subroutine Cool_init
