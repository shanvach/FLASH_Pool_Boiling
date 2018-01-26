module Plasma_interface

        implicit none

#include "constants.h"
#include "Flash.h"

        interface
        subroutine Plasma_init(blockCount,blockList,restart)
                implicit none
                integer, INTENT(INOUT) :: blockCount
                integer, INTENT(INOUT) :: blockList(MAXBLOCKS)
                logical, INTENT(IN)    :: restart
        end subroutine Plasma_init
        end interface

        interface
        subroutine Plasma_finalize()
                implicit none
        end subroutine
        end interface

        interface
        subroutine Plasma(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)
                implicit none
                integer, INTENT(INOUT) :: blockCount
                integer, INTENT(IN) :: sweepOrder
                integer, INTENT(INOUT) :: blockList(MAXBLOCKS)
                real,    INTENT(IN) :: timeEndAdv, dt, dtOld
        end subroutine Plasma
        end interface


        interface
        subroutine Plasma_Solve(T_p, T_gen, T_o, dfun, dcoeff, dt, dx, dy, ix1,ix2, jy1, jy2,T_res)
                implicit none
                real, dimension(:,:,:), intent(inout) :: T_p
                real, dimension(:,:,:), intent(in) :: T_o, T_gen, dfun, dcoeff
                real, intent(in) :: dt, dx, dy
                integer, intent(in) :: ix1, ix2, jy1, jy2
                real, intent(out) :: T_res
        end subroutine
        end interface

        interface
        subroutine Plasma_get_reactions(Te,Th)
                implicit none
                real, intent(in) :: Te, Th
        end subroutine Plasma_get_reactions
        end interface

        interface
        subroutine Plasma_computeDt(pls_mindt,pls_minloc)
                implicit none
                real, intent(INOUT) :: pls_mindt
                integer, intent(INOUT) :: pls_minloc(5)
        end subroutine Plasma_computeDt
        end interface

        interface
        subroutine Plasma_computeDtLocal(blockID,   &
                              isize, jsize, ksize,  &
                              dx, dy, dz,           &
                              blkLimits,blkLimitsGC,&
                              solnData,&
                              facexData,faceyData,  &
                              facezData,            &
                              dtLocal, lminloc)
                implicit none
                integer, intent(IN) :: blockID
                integer,dimension(2,MDIM), intent(IN) :: blkLimits,blkLimitsGC
                integer, intent(IN) :: isize,jsize,ksize
                real, intent(IN) :: dx, dy, dz
                real, pointer,dimension(:,:,:,:)  :: solnData,facexData,faceyData,facezData
                real, intent(INOUT) :: dtLocal
                integer, intent(INOUT) :: lminloc(5)

        end subroutine Plasma_computeDtLocal
        end interface
        
        interface
        subroutine Plasma_hvDiffCoeff( DiffCoeff, P_h, T_h, ix1, ix2, jy1, jy2,RSCD, MHSP)
                implicit none
                real, dimension(:,:,:), intent(inout) :: DiffCoeff
                real, dimension(:,:,:), intent(in) :: P_h, T_h
                real, intent(in) :: RSCD, MHSP
                integer, intent(in) :: ix1, ix2, jy1, jy2       
        end subroutine Plasma_hvDiffCoeff 
        end interface

        interface
        subroutine Plasma_elDiffCoeff(Diffcoeffion, DiffCoeffel, T_e, T_h,&
                                      Fvea, Fvei, ix1, ix2,jy1, jy2)
                implicit none
                real, dimension(:,:,:), intent(inout) :: DiffCoeffion,DiffCoeffel
                real, dimension(:,:,:), intent(in) :: T_e, T_h, Fvea, Fvei
                integer, intent(in) :: ix1, ix2, jy1, jy2
        end subroutine Plasma_elDiffCoeff
        end interface

        interface
        subroutine Plasma_ColFreq(vei, vea, N_a, N_e, T_e, ix1, ix2, jy1, jy2)
                implicit none
                real, dimension(:,:,:), intent(in) :: N_a, N_e, T_e
                real, dimension(:,:,:), intent(inout) :: vea, vei
                real :: Kel, Lna
                integer, intent(in) :: ix1, ix2, jy1, jy2
        end subroutine Plasma_ColFreq
        end interface
        
        interface
        subroutine Plasma_sumNeutrals(N_as, N_at, ix1, ix2, jy1, jy2)
                implicit none
                real, dimension(:,:,:), intent(in) :: N_as
                real, dimension(:,:,:), intent(inout) :: N_at
                integer, intent(in) :: ix1, ix2, jy1, jy2
        end subroutine Plasma_sumNeutrals
        end interface

        interface
        subroutine Plasma_sumIons(N_is, N_it, ix1, ix2, jy1, jy2 )
                implicit none
                real, dimension(:,:,:), intent(in) :: N_is
                real, dimension(:,:,:), intent(inout) :: N_it
                integer, intent(in) :: ix1, ix2, jy1, jy2
        end subroutine Plasma_sumIons
        end interface

        interface
        subroutine Plasma_spReactions(RSP0,RSP1,RSP2,RSP3,RSP4,RSP5,RSP6,RSP7,&
                                      RSP8,RSP9,RSP10,RSP11,RSP12,RSP13,RSP14,& 
                                      RSP15,RSP16,T_h,T_e,ix1, ix2, jy1, jy2)
                implicit none
                real, dimension(:,:,:), intent(inout) :: RSP0,RSP1,RSP2,RSP3,RSP4,&
                                                         RSP5,RSP6,RSP7,RSP8,RSP9,&
                                                         RSP10,RSP11,RSP12,RSP13,&
                                                         RSP14,RSP15,RSP16
                real, dimension(:,:,:), intent(in) :: T_h, T_e
                integer, intent(in) :: ix1, ix2, jy1, jy2
        end subroutine Plasma_spReactions
        end interface
   
        interface
        subroutine Plasma_spGeneration(N_h0,N_h1,N_h2,N_h3,N_h4,N_h5,N_h6,& 
                                       N_h7,N_h8,N_h9,N_el,RSP0,RSP1,RSP2, &
                                       RSP3,RSP4,RSP5,RSP6,RSP7,RSP8,RSP9,&
                                       RSP10,RSP11,RSP12,RSP13,RSP14,RSP15,&
                                       RSP16,GNH0,GNH1,GNH2,GNH3,GNH4,GNH5,&
                                       GNH6,GNH7,GNH8,GNH9,GNE,GNEBZ,GNERT,&
                                       ix1,ix2,jy1,jy2)
                implicit none
                real, dimension(:,:,:), intent(in) :: N_h0,N_h1,N_h2,N_h3,N_h4,&
                                                      N_h5,N_h6,N_h7,N_h8,N_h9,&
                                                      N_el
                real, dimension(:,:,:), intent(in) :: RSP0,RSP1,RSP2,RSP3,RSP4,&
                                                      RSP5,RSP6,RSP7,RSP8,RSP9,&
                                                      RSP10,RSP11,RSP12,RSP13,&
                                                      RSP14,RSP15,RSP16
                real, dimension(:,:,:), intent(inout) :: GNH0,GNH1,GNH2,GNH3,&
                                            GNH4,GNH5,GNH6,GNH7,&
                                            GNH8,GNH9,GNE,GNEBZ,&
                                            GNERT
                integer, intent(in) :: ix1, ix2, jy1, jy2
        end subroutine Plasma_spGeneration
        end interface

end module Plasma_interface
