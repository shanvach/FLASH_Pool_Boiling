!!****if* source/Simulation/SimulationMain/AdvectMassScalars/wr_integrals
!!
!!  NAME
!!    wr_integrals
!!
!!  SYNOPSIS
!!    call wr_integrals(isfirst, simtime)
!!    call wr_integrals(integer, real)
!!
!!  DESCRIPTION
!!     Compute the values of integral quantities (eg. total energy)
!!     and write them to an ASCII file.  If this is the initial step,
!!     create the file and write a header to it before writing the
!!     data.
!!
!!     Presently, this supports 1, 2, and 3-d Cartesian geometry and 2-d
!!     cylindrical geometry (r,z).  More geometries can be added by modifying
!!     the volume of each zone (dvol) computed in SumLocalIntegrals.
!!
!!***

      subroutine WriteIntHeader (funit)

!=============================================================================

        use multifluid_database, ONLY : n_fluids, get_fluid_property

implicit none !! Added by fix script
        character(len=80) :: fluid_name(n_fluids)

        integer :: funit, i

!=============================================================================

        do i = 1, n_fluids
           call get_fluid_property( i, "name", fluid_name(i) )
        end do

        write (funit, 10)               &
           '#time                    ', &
           'mass                     ', &
           'x-momentum               ', &
           'y-momentum               ', & 
           'z-momentum               ', &
           'E_total                  ', &
           'E_internal               ', &
           'E_kinetic                ', &
           'mass (ms_1)              ', &
           'mass (ms_2)              ', &
           'mass (ms_3)              ', &
           'mass (ms_4)              ', &
           'mass (ms_5)              ', &
           (fluid_name(i)(1:min(len(fluid_name(i)),22)),i=1,n_fluids)
      
10      format (2x,50(a22, :, 1X))

!=============================================================================

        return
      end subroutine WriteIntHeader
