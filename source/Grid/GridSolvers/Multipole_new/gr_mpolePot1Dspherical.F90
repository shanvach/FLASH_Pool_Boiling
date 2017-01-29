!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpolePot1Dspherical
!!
!! NAME
!!
!!  gr_mpolePot1Dspherical
!!
!! SYNOPSIS
!!
!!  gr_mpolePot1Dspherical  (integer, intent(in) :: ipotvar)
!!
!! DESCRIPTION
!!
!!  Computes the potential field for a one-dimensional spherical geometry
!!  using the mass moments already calculated. On output the variable
!!  indexed by ipotvar contains the potential. The calculations are
!!  entirely local to each processor, since each processor has a local
!!  copy of the moments.
!!
!! ARGUMENTS
!!
!!  ipotvar : index to variable containing the potential
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpolePot1Dspherical (ipotvar)

  use Grid_interface,    ONLY : Grid_getBlkPtr,        &
                                Grid_releaseBlkPtr,    &
                                Grid_getBlkBoundBox,   &
                                Grid_getDeltas,        &
                                Grid_getBlkIndexLimits

  use gr_mpoleData,      ONLY : gr_mpoleGravityConstant,        &
                                gr_mpoleDrInv,                  &
                                gr_mpoleDrInnerZoneInv,         &
                                gr_mpoleMaxRadialZones,         & 
                                gr_mpoleMinRadialZone,          & 
                                gr_mpoleZoneRmax,               &
                                gr_mpoleZoneQmax,               &
                                gr_mpoleZoneType,               &
                                gr_mpoleZoneScalarInv,          &
                                gr_mpoleZoneLogNormInv,         &
                                gr_mpoleZoneExponentInv,        &
                                gr_mpoleInnerZoneMaxR,          &
                                gr_mpoleInnerZoneDrRadii,       &
                                gr_mpoleInnerZoneQlower,        &
                                gr_mpoleInnerZoneQupper,        &
                                gr_mpoleInnerZoneResolution,    &
                                gr_mpoleInnerZoneResolutionInv, &
                                gr_mpoleOuterZoneQshift,        &
                                gr_mpoleQDampingI,              &
                                gr_mpoleMomentR,                &
                                gr_mpoleMomentI,                &
                                gr_mpoleBlockCount,             &
                                gr_mpoleBlockList

  implicit none
  
#include "Flash.h"
#include "constants.h"
#include "gr_mpole.h"
  
  integer, intent (in) :: ipotvar

  logical :: i2
  logical :: innerZonePotential

  integer :: blockNr, blockID
  integer :: DrUnit
  integer :: imax, imin, iC, iCmax, iF, iFmax
  integer :: Q, Qlocal, Qlower, Qupper
  integer :: type
  integer :: zone

  integer :: blkLimits   (LOW:HIGH,1:MDIM)
  integer :: blkLimitsGC (LOW:HIGH,1:MDIM)

  real    :: bndBoxILow
  real    :: DeltaI, DeltaIHalf
  real    :: facePotential
  real    :: Qfloat, QfracI, QfracR
  real    :: RdotI, IdotR
  real    :: rlocal, rinDrs, rinvI
  real    :: Rsph
  real    :: sclInv, lgnInv, expInv

  real    :: delta           (1:MDIM)
  real    :: bndBox (LOW:HIGH,1:MDIM)

  real, pointer :: solnData (:,:,:,:)
!  
!
!     ...Sum quantities over all locally held leaf blocks.
!
!
!$omp do schedule (static)
  do blockNr = 1,gr_mpoleBlockCount

     blockID = gr_mpoleBlockList (blockNr)

     call Grid_getBlkBoundBox     (blockID, bndBox)
     call Grid_getDeltas          (blockID, delta)
     call Grid_getBlkPtr          (blockID, solnData)
     call Grid_getBlkIndexLimits  (blockID, blkLimits, blkLimitsGC)

     imin       = blkLimits (LOW, IAXIS)
     imax       = blkLimits (HIGH,IAXIS)

     iCmax      = imax - imin + 1      ! # of local (inner) cells in i direction
     iFmax      = iCmax                !  max local (inner) face index in i direction (1st face -> index 0)

     DeltaI     = delta (IAXIS)
     DeltaIHalf = DeltaI * HALF

     bndBoxILow = bndBox (LOW,IAXIS)

     solnData (ipotvar , imin:imax , 1,1) = ZERO
!
!
!          ...The 1D spherical case. In this case each point is characterized by
!             its radius Rsph from the center (stored in the i-index). Only L = 0
!             solid harmonics are calculated and combined with the moments.
!
!
     Rsph = bndBoxILow

     do iF = 0, iFmax                                   ! loop over local i face indices

      iC = iF + 1                                       ! local (inner) largest i cell index for i face
      i2 = iC > 1 .and. iC < iCmax + 1                  ! 2 cells i and i-1 share i face?
      iC = imin - 1 + min (iC , iCmax)                  ! change to global (inner + guard) i cell index
!
!
!        ...Find the radial bin.
!
!
      innerZonePotential = Rsph <= gr_mpoleInnerZoneMaxR


      if (innerZonePotential) then

          rinDrs = Rsph * gr_mpoleDrInnerZoneInv
          DrUnit = int (ceiling (rinDrs))
          Qlower = gr_mpoleInnerZoneQlower (DrUnit)
          Qupper = gr_mpoleInnerZoneQupper (DrUnit)
          QfracR = ZERO
          QfracI = ONE

          do Q = Qlower,Qupper
             if (rinDrs <= gr_mpoleInnerZoneDrRadii (Q)) exit
          end do

      else

          do zone = gr_mpoleMinRadialZone , gr_mpoleMaxRadialZones
             if (Rsph - gr_mpoleZoneRmax (zone) <= ZERO) exit
          end do

          rlocal = Rsph - gr_mpoleZoneRmax (zone - 1)
          type   = gr_mpoleZoneType        (zone)
          sclInv = gr_mpoleZoneScalarInv   (zone)
          expInv = gr_mpoleZoneExponentInv (zone)

          if (type == ZONE_EXPONENTIAL) then
              Qfloat = (rlocal * sclInv * gr_mpoleDrInv) ** expInv
          else if (type == ZONE_LOGARITHMIC) then
              lgnInv = gr_mpoleZoneLogNormInv (zone)
              Qfloat = expInv * log (rlocal * sclInv * gr_mpoleDrInv * lgnInv + ONE)
          end if

          Qlocal = ceiling (Qfloat)
          QfracI = real (Qlocal) - Qfloat
          QfracR = ONE - QfracI
          Q      = gr_mpoleZoneQmax (zone - 1) + Qlocal + gr_mpoleOuterZoneQshift

      end if
!
!
!        ...Calculate and add the current face potential to the relevant cell(s) of the
!           potential block. Note, that a face can only be shared by up to 2 cells and
!           no more.
!
!
      rinvI = ONE / Rsph

      RdotI =   QfracI * gr_mpoleQDampingI (Q)   * gr_mpoleMomentI (1,Q  )   &
              + QfracR * gr_mpoleQDampingI (Q+1) * gr_mpoleMomentI (1,Q+1)

      IdotR = rinvI * (  QfracI * gr_mpoleMomentR (1,Q-1) &
                       + QfracR * gr_mpoleMomentR (1,Q)   )

      facePotential = - gr_mpoleGravityConstant * (RdotI + IdotR)

      if (i2) then
          solnData (ipotvar,iC-1,1,1) = solnData (ipotvar,iC-1,1,1) + facePotential
          solnData (ipotvar,iC  ,1,1) = solnData (ipotvar,iC  ,1,1) + facePotential
      else
          solnData (ipotvar,iC,1,1) = solnData (ipotvar,iC,1,1) + facePotential
      end if

      Rsph = Rsph + DeltaI

     end do
!
!
!    ...Form the potential average in each cell.
!
!
     solnData (ipotvar,imin:imax,1,1) = HALF * solnData (ipotvar,imin:imax,1,1)
!
!
!    ...Get ready for retrieving next LEAF block for the current processor.
!
!
     call Grid_releaseBlkPtr (blockID, solnData)

  end do
!$omp end do
!
!
!    ...Ready!
!
!
  return
end subroutine gr_mpolePot1Dspherical

