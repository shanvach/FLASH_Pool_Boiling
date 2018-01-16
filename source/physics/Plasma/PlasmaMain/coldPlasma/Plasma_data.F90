module Plasma_data

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


end module Plasma_data
