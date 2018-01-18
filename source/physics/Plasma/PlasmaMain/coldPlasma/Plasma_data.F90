module Plasma_data

! declare variables, which can be used as constants during sim.
! for example: electron mass, charge, ion masses  which do not change with time.

        implicit none
        
        integer,save :: pls_cflflg
        logical,save :: pls_restart

        integer, save :: pls_nstep
        integer, save :: pls_meshMe
        integer, save :: pls_meshNumProcs
        integer, save :: pls_meshComm

        real, save    :: pls_cfl
        real, save    :: pls_sigma
        real, save    :: pls_dtspec


        logical, save :: pls_predcorrflg

        logical, save :: pls_outflowgridChanged

        integer, save :: pls_prol_method

        real, save :: pls_dcoeff

        real, dimension(10), save :: pls_RSCD ! Collision Diameter, Angstrom
        real, dimension(10), save :: pls_MHSP ! binary molar mass g/mol

        real, dimension(10), save :: pls_MMIX ! molar mass
        
        real, save :: pls_Ckb ! Boltzmann Constant J/K
        real, save :: pls_Cme ! electron mass kg
        real, save :: pls_Ce  ! elementary charge, C
        real, save :: pls_Cpi ! pi
        real, save :: pls_gam ! Euler's number

end module Plasma_data
