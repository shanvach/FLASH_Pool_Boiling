!!****if* source/physics/materialProperties/Opacity/OpacityMain/save/Opacity_unitTest
!!
!! NAME
!!
!!  Opacity_unitTest
!!
!! SYNOPSIS
!!
!!  Opacity_unitTest (integer (in)    :: fileUnit,
!!                    logical (inout) :: perfect)
!!
!! DESCRIPTION
!!
!!  This is a unitTest setup for testing the Opacity unit. Normally called
!!  from Driver_evolveFlash within a Simulation unitTest. See for example
!!  the directory 'source/Simulation/SimulationMain/unitTest/Opacity'.
!!
!!  Performs a test run on the routines that interpolate the generated opacity
!!  tables. Synthetic opacity tables have been produced during the unitTest
!!  initialization according to the number of species.
!!
!! ARGUMENTS
!!
!!  fileUnit : number of file unit for diagnostic output
!!  perfect  : will be set .true., if the test is correct
!! 
!! NOTES
!!
!!  In the Simulation unit, you must set your Config file to 
!!  REQUIRES Opacity/OpacityMain/unitTest
!!
!!
!!***
subroutine Opacity_unitTest (fileUnit, perfect)

  use Opacity_dataUnitTest, ONLY : op_totalSpecies,          &
                                   op_absorptionKind,        &
                                   op_emissionKind,          &
                                   op_transportKind,         &
                                   op_nstepsTemperature,     &
                                   op_nstepsDensity,         &
                                   op_temperatureFirst,      &
                                   op_temperatureStep,       &
                                   op_ionNumberDensityFirst, &
                                   op_ionNumberDensityStep,  &
                                   op_log10opacityFirstPA,   &
                                   op_log10opacityFirstPE,   &
                                   op_log10opacityFirstRO,   &
                                   op_log10opacityStepPA,    &
                                   op_log10opacityStepPE,    &
                                   op_log10opacityStepRO,    &
                                   op_speciesWeights,        &
                                   op_massFractions,         &
                                   op_Avogadro,              &
                                   zero,one,two,ten

  use Driver_interface,     ONLY : Driver_abortFlash
  use Opacity_interface,    ONLY : Opacity
  use Grid_interface,       ONLY : Grid_getListOfBlocks,    &
                                   Grid_getBlkPtr,          &
                                   Grid_releaseBlkPtr,      &
                                   Grid_getBlkIndexLimits

  use op_interface,         ONLY : op_LaguerreQuadratureRule

  implicit none

# include "constants.h"
# include "Flash.h"

  integer, intent (in)    :: fileUnit
  logical, intent (inout) :: perfect

  character (len=20) :: kindA,kindE,kindT

  integer :: blkCount
  integer :: block
  integer :: blockID
  integer :: cell
  integer :: cellIndex
  integer :: firstCellIndex
  integer :: species
  integer :: totalCells

  real    :: speciesWeight
  real    :: cellmassDensity
  real    :: cellTemperature
  real    :: cT,cD
  real    :: log10opacityFirstPA
  real    :: log10opacityFirstPE
  real    :: log10opacityFirstRO
  real    :: log10opacityStepPA
  real    :: log10opacityStepPE
  real    :: log10opacityStepRO
  real    :: massFraction
  real    :: nT,nD
  real    :: OA,OE,OT
  real    :: opacityAbsorption
  real    :: opacityEmission
  real    :: opacityTransport
  real    :: sumOA,sumOE,sumOT,sumMW
  real    :: Tf,Ts,Df,Ds

  integer, dimension (MAXBLOCKS) :: blkList

  integer, dimension (2,MDIM)    :: blkLimits
  integer, dimension (2,MDIM)    :: blkLimitsGC

  real, pointer, dimension (:,:,:,:) :: solnData




  real    :: Integral
  real    :: r,w
  real    :: T
  real    :: beta
  integer :: i,n
  integer :: nRoots
  real, allocatable :: Roots   (:)
  real, allocatable :: Weights (:)

!
!
!   ...Test the Biggs stuff.
!
!  
!  call op_readBiggs1971xraydatFile ()
!
!  return
!

  call op_initLowTemp ()
  call op_initIntegrate ()
  call op_PlanckMeanOpacity ()
  call op_finalizeIntegrate ()

  return

!  call op_setAtomNames ()
!  call op_setPEarrayJmax ()
!  call op_setPEenergyRange ()
!  call op_setPEcoeffsAij4 ()
!  call op_writePEdata ()
!  call op_writeAtomPEopacity2file (55,1.0,100.0,1000)
!  call op_finalizeLowTemp ()
!
!
!   ...Test the Quadrature stuff.
!
!  
   write (*,*) ' beta ? '
   read  (*,*) beta
   write (*,*) ' Integration limit T ? '
   read  (*,*) T
   write (*,*) ' # of Roots ? '
   read  (*,*) nRoots

   allocate (Roots   (1:nRoots))
   allocate (Weights (1:nRoots))

   call op_initIntegrate ()


   do n = 1,nRoots
      call op_LaguerreQuadratureRule (n,beta-1.0,T,Roots,Weights)

      Integral = 0.0
      do i = 1,n
         r = Roots (i)
         w = Weights (i)
         Integral = Integral + (r*exp (r)/(exp(r)-1.0)) * w
      end do

      write (*,*) ' n,Integral = ',n,Integral
   end do

   call op_finalizeIntegrate ()

   deallocate (Roots)
   deallocate (Weights)

   return
!
!
!   ...Get all leaf blocks.
!
!  
  call Grid_getListOfBlocks (LEAF, blkList, blkCount)

  do block = 1,blkCount

     blockID = blkList (block)
     call Grid_getBlkPtr         (blockID, solnData)
     call Grid_getBlkIndexLimits (blockID, blkLimits, blkLimitsGC)

     firstCellIndex = blkLimits (LOW,IAXIS)

     totalCells = blkLimits (HIGH,IAXIS) - blkLimits (LOW,IAXIS) + 1

     do cell = 1,totalCells

        cellIndex       = firstCellIndex + cell - 1
        cellmassDensity = solnData (DENS_VAR,cellIndex,1,1)
        cellTemperature = solnData (TELE_VAR,cellIndex,1,1)

        write (*,*) ' cell # ',cell
        write (*,*) ' ---------------------------------------- '
        write (*,*) ' massDensity in cell          = ',cellmassDensity
        write (*,*) ' electronTemperature in cell  = ',cellTemperature
!
!
!   ...Get the tabulated opacities.
!
!  
        call Opacity (solnData (1:NUNK_VARS,cellIndex,1,1), &
                      1,                          &
                      opacityAbsorption,          &
                      opacityEmission,            &
                      opacityTransport            )

        write (*,*) ' opacityAbsorption (tabulated) = ',opacityAbsorption
        write (*,*) ' opacityEmission   (tabulated) = ',opacityEmission
        write (*,*) ' opacityTransport  (tabulated) = ',opacityTransport
!
!
!   ...Calculate the analytic opacities.
!
!  
        sumOA = zero
        sumOE = zero
        sumOT = zero
        sumMW = zero

        do species = 1,op_totalSpecies

           kindA = op_absorptionKind (species)
           kindE = op_emissionKind   (species)
           kindT = op_transportKind  (species)

           speciesWeight = op_speciesWeights (species)
           massFraction  = op_massFractions  (species)

           cT = log10 (cellTemperature)
           cD = log10 (massFraction * cellmassDensity * op_Avogadro / speciesWeight)
           Tf = log10 (op_temperatureFirst      (species))
           Ts = log10 (op_temperatureStep       (species))
           Df = log10 (op_ionNumberDensityFirst (species))
           Ds = log10 (op_ionNumberDensityStep  (species))

           nT = one + (cT - Tf) / Ts
           nD = one + (cD - Df) / Ds

           log10opacityFirstPA = op_log10opacityFirstPA (species)
           log10opacityFirstPE = op_log10opacityFirstPE (species)
           log10opacityFirstRO = op_log10opacityFirstRO (species)
           log10opacityStepPA  = op_log10opacityStepPA  (species)
           log10opacityStepPE  = op_log10opacityStepPE  (species)
           log10opacityStepRO  = op_log10opacityStepRO  (species)
!
!
!   ...Calculate the absorption opacity for the current species.
!
!  
           if (kindA == "Planck Absorption") then
               OA = ten ** (log10opacityFirstPA + (nT + nD - two) * log10opacityStepPA)
           else if (kindA == "Planck Emission") then
               OA = ten ** (log10opacityFirstPE + (nT + nD - two) * log10opacityStepPE)
           else if (kindA == "Rosseland") then
               OA = ten ** (log10opacityFirstRO + (nT + nD - two) * log10opacityStepRO)
           else
               call Driver_abortFlash ('[Opacity_unitTest] ERROR: Absorption opacity kind not recognized')
           end if
!
!
!   ...Calculate the emission opacity for the current species.
!
!  
           if (kindE == "Planck Absorption") then
               OE = ten ** (log10opacityFirstPA + (nT + nD - two) * log10opacityStepPA)
           else if (kindE == "Planck Emission") then
               OE = ten ** (log10opacityFirstPE + (nT + nD - two) * log10opacityStepPE)
           else if (kindE == "Rosseland") then
               OE = ten ** (log10opacityFirstRO + (nT + nD - two) * log10opacityStepRO)
           else
               call Driver_abortFlash ('[Opacity_unitTest] ERROR: Emission opacity kind not recognized')
           end if
!
!
!   ...Calculate the transport opacity for the current species.
!
!  
           if (kindT == "Planck Absorption") then
               OT = ten ** (log10opacityFirstPA + (nT + nD - two) * log10opacityStepPA)
           else if (kindT == "Planck Emission") then
               OT = ten ** (log10opacityFirstPE + (nT + nD - two) * log10opacityStepPE)
           else if (kindT == "Rosseland") then
               OT = ten ** (log10opacityFirstRO + (nT + nD - two) * log10opacityStepRO)
           else
               call Driver_abortFlash ('[Opacity_unitTest] ERROR: Transport opacity kind not recognized')
           end if
!
!
!   ...Calculate .
!
!  
           sumOA = sumOA + (massFraction * OA) / speciesWeight
           sumOE = sumOE + (massFraction * OE) / speciesWeight
           sumOT = sumOT + (massFraction * OT) / speciesWeight
           sumMW = sumMW +  massFraction       / speciesWeight

           write (*,*) ' speciesWeight = ',speciesWeight
           write (*,*) ' massFraction  = ',massFraction
           write (*,*) ' OA = ',OA
           write (*,*) ' OE = ',OE
           write (*,*) ' OT = ',OT


        end do

        opacityAbsorption = cellmassDensity * sumOA / sumMW
        opacityEmission   = cellmassDensity * sumOE / sumMW
        opacityTransport  = cellmassDensity * sumOT / sumMW

        write (*,*) ' opacityAbsorption  (analytic) = ',opacityAbsorption
        write (*,*) ' opacityEmission    (analytic) = ',opacityEmission
        write (*,*) ' opacityTransport   (analytic) = ',opacityTransport

     end do

     call Grid_releaseBlkPtr (blockID, solnData)

  end do

  perfect = .true.
!
!
!    ...Ready!
!
!
  return
end subroutine Opacity_unitTest
