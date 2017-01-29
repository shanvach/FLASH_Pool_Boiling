!!****if* source/Simulation/SimulationMain/Nuc2Grid/Grid_convolve
!!
!! NAME
!!
!!  Grid_convolve
!!
!!
!! SYNOPSIS
!!
!!  call Grid_convolve(integer(IN) :: iSoln,
!!                     integer(IN) :: iSrc, 
!!                     integer(IN) :: iKern, 
!!                  integer(6)(IN) :: bcTypes,
!!                   real(2,6)(IN) :: bcValues,
!!                     real(INOUT) :: scaleFactor)
!!
!! DESCRIPTION
!!
!!   Convolution routine.  This interface module provides the 
!!   fft based method of convolution.
!!
!!   Initially implemented for periodic boundary conditions
!!   on a uniform grid.  
!!
!!
!! ARGUMENTS
!!
!!  iSoln -  index to variable containing potential (smeared distribution)
!!  iSrc  - index to variable containing input (distribution to be smeared)
!!  iKern - index to variable containing the convolution Kernel
!!  bcTypes - boundary types along various faces,
!!          valid values are: (although only some are implemented)
!!          PERIODIC -- supported
!!          DIRICHLET
!!          OUTFLOW (Neuman)
!!          HYDROSTATIC (Dirichlet value)
!!          ISOLATED
!!  bcValues - the values to boundary conditions, currently not used
!!  scaleFactor -  factor used to scale results
!!
!!***

subroutine Grid_convolve (iSoln, iSrc, iKern, bcTypes, bcValues, scaleFactor)

  use gr_pfftData, ONLY : pfft_inLen,pfft_outLen,pfft_setupOnce,pfft_usableProc

  use Grid_interface, ONLY : Grid_pfftMapToInput,Grid_getGlobalIndexLimits,&
       Grid_getBlkIndexLimits,Grid_pfftInit, Grid_pfft,Grid_pfftMapFromOutput, &
       Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_pfftFinalize
  use gr_pfftInterface, ONLY : gr_pfftSpecifyTransform
  use Grid_data, ONLY : gr_meshNumProcs, gr_meshMe
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Pfft.h"  

  integer, intent(in)    :: iSoln, iSrc, iKern
  integer, intent(in)    :: bcTypes(2*MDIM)
  real, intent(in)       :: bcValues(2,2*MDIM)
  real, intent(inout)    :: scaleFactor

  !--------------------------------------------------------------------------
  real, allocatable, dimension(:) :: inArray,tranArray,outArray
  real, allocatable, dimension(:),SAVE :: cKernArray
  real :: cKernScale

  integer, dimension(MDIM) :: localSize,globalSize,transformType
  integer, dimension(0:MDIM) :: baseDatType
  integer :: inSize,tranSize,i
  logical :: needMap
  logical,save :: kernelTransformed = .FALSE.


  if(.not.pfft_setupOnce) then
     needMap= (NDIM > 1)   ! For 1D, PFFT requires to be run on only 1 proc anyway.
     call Grid_getGlobalIndexLimits(globalSize)

     call gr_pfftSpecifyTransform(transformType, baseDatType)

     call Grid_pfftInit(NDIM,needMap,globalSize,&
          localSize,transformType, baseDatType)

     
  end if

  !Important.  Tests that this processor should be doing work
  if(.not.pfft_usableProc) return   

  inSize=pfft_inLen(IAXIS)*pfft_inLen(JAXIS)*pfft_inLen(KAXIS)
  tranSize=2*pfft_outLen(IAXIS)*pfft_outLen(JAXIS)*pfft_outLen(KAXIS)

  allocate(inArray(inSize+2))

  if (.not.kernelTransformed) then
     allocate(cKernArray(transize))
     print*,'transize=',transize
     cKernArray(:) = -1.0
     call Grid_pfftMapToInput(iKern,inArray)  
     kernelTransformed = .TRUE.
     ! Forward transform of density 
     call Grid_pfft(PFFT_FORWARD,inArray,cKernArray)
     cKernScale = pfft_outLen(IAXIS)*pfft_outLen(JAXIS)*pfft_outLen(KAXIS)
     cKernArray(:) = cKernArray(:) * cKernScale
     do i=1,transize/2-1,2
999     format(I4,':',2(1PG20.12))
!        print 999,i,cKernArray(i:i+1)
     end do
  end if

  allocate(outArray(inSize+2))
  allocate(tranArray(tranSize))


  !! Here's the real work of the fft
  ! Converts to uniform mesh (on output, inArray contains uniformly mapped density)
  call Grid_pfftMapToInput(iSrc,inArray)  
  ! Forward transform of density 
  call Grid_pfft(PFFT_FORWARD,inArray,tranArray)
  ! Calculates the transform of iSoln = GPOT
  !  which is the transform of the delSquared(u) = rho in Poisson equation
!!$  if (iSoln==NUMP_VAR) then
!!$     print*
!!$     do i=1,transize-1,2
!!$        print 999,i,tranArray(i:i+1)
!!$     end do
!!$  end if
  call gr_pfftMult(tranArray, cKernArray, .FALSE.)
!!$  if (iSoln==NUMP_VAR) then
!!$     print*
!!$     do i=1,transize-1,2
!!$        print 999,i,tranArray(i:i+1)
!!$     end do
!!$  end if
  ! Inverse transform of GPOT
  call Grid_pfft(PFFT_INVERSE,tranArray,outArray)
  ! Now multiply by the poisson factor
   outArray(1:inSize) = outArray(1:inSize)*scaleFactor

  ! Map back to the non-uniform mesh
  call Grid_pfftMapFromOutput(iSoln,outArray)


  deallocate(inArray)
  deallocate(tranArray)
  deallocate(outArray)
  if(.not.pfft_setupOnce) then
     call Grid_pfftFinalize()
  end if

  return
end subroutine Grid_convolve
