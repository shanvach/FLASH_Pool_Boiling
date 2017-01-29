!******************************************************************************

!  Routine:     init_materials()

!  Description: Initializes the materials module, in particular the multifluid
!               database and the equation of state.  This version sets up
!               three tracer fluids for the three-layer target problem.


        subroutine init_materials

!==============================================================================

        use dBase, ONLY : dbasePropertyInteger
        use multifluid_database

        implicit none

        integer urp
        integer num_ions
!==============================================================================

        num_ions = dBasePropertyInteger("NumberOfSpecies")

!               Set up the multifluid database.

        call init_mfluid_db (num_ions)          
                                                

        call add_fluid_to_db ("cu ", "cu ", A=63.54e0, & 
     &                         Z=29.0e0, gamma=2.0e0)
        call add_fluid_to_db ("ch ", "ch ", A=13.e0, & 
     &                         Z=7.0e0, gamma=2.0e0)
        call add_fluid_to_db ("cf", "cf", A=12.0e0, & 
     &                          Z=6.0e0, gamma=1.3e0)


        call find_fluid_index ("cu ", urp)
        write (*,*) "cu  = ", urp
        call find_fluid_index ("ch", urp)
        write (*,*) "ch = ", urp
        call find_fluid_index ("cf", urp)
        write (*,*) "cf = ", urp

!==============================================================================

        return
        end

