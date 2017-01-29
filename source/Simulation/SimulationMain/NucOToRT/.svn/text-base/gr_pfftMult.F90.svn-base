!!****if* source/Simulation/SimulationMain/NucOToRT/gr_pfftMult
!!
!! NAME 
!!
!!   gr_pfftMult
!!
!! SYNOPSIS
!!
!!   call gr_pfftMult(real(INOUT) :: transArray(*), 
!!                    real(IN)    :: factorArray(*),
!!                    logical(IN) :: conjugate)
!!
!! DESCRIPTION 
!!
!!  Multiply two functions in momentum space.
!! 
!! ARGUMENTS
!!
!!  transArray : At input it contains transformed data, at 
!!               output it has the product of the transformed data with the
!!               factor array. 
!!  factorArray : At input it contains a factor array, which is also in
!!                transformed space.
!!  conjugate :  flag that indicates whether the complex conjugate of factorArray
!!               should be used,
!!               rather than using factorArray directly.
!!
!!
!!***
subroutine gr_pfftMult(transArray, factorArray, conjugate)

#include "constants.h"  
#include "Pfft.h"

  use gr_pfftData, ONLY : pfft_wave, pfft_localLimits, pfft_outLen,&
       pfft_transformType,pfft_dimOrder,pfft_ndim, pfft_usableProc

  use Grid_interface, ONLY : Grid_getDeltas

  implicit none
  
  real,dimension(*),intent(INOUT) :: transArray
  real,dimension(*),intent(IN) :: factorArray
  logical :: conjugate

  integer :: jfactor, kfactor
  real,dimension(MDIM) :: delta
  integer :: blockID=1
  real :: partialJ,partialK,partial
  integer :: i,j,k,ii,jj,kk,n
  logical :: isComplex, isPacked, packedReals
  integer :: multiplier,jOffset,kOffset

  if (.not.pfft_usableProc) return

  isComplex=pfft_transformType(pfft_dimOrder(pfft_ndim))==PFFT_COMPLEX
  isComplex=(pfft_transformType(pfft_dimOrder(pfft_ndim))==PFFT_REAL2C)&
       .or.isComplex
  isComplex=(pfft_transformType(pfft_dimOrder(pfft_ndim))==PFFT_REAL2C_STUFF)&
       .or.isComplex
  isComplex=(pfft_transformType(pfft_dimOrder(pfft_ndim))==PFFT_REAL2C_EXTEND)&
       .or.isComplex
  if (isComplex) then
     isPacked = (pfft_transformType(1)==PFFT_REAL2C_STUFF)
     if(.NOT.isPacked .AND. pfft_ndim > 1)&
          isPacked = (pfft_transformType(2)==PFFT_REAL2C_STUFF)
     if(.NOT.isPacked .AND. pfft_ndim > 2)&
          isPacked = (pfft_transformType(3)==PFFT_REAL2C_STUFF)
  else
     isPacked = .FALSE.
  end if

  multiplier=1
  if(isComplex)then
     multiplier=2
  end if

  jfactor=max(1,pfft_outLen(IAXIS))*multiplier
  kfactor=max(1,jfactor*pfft_outLen(JAXIS))

  kk=0
  kOffset=0
  do k = pfft_localLimits(LOW,KAXIS),pfft_localLimits(HIGH,KAXIS)
     if(pfft_ndim==MDIM) then
        kk=kk+1
        partialK=pfft_wave(kk,KAXIS)
     else
        partialK=0.0
     end if
     jj=0
     jOffset=0
     do j = pfft_localLimits(LOW,JAXIS),pfft_localLimits(HIGH,JAXIS)
        if(pfft_ndim>1) then
           jj=jj+1
           partialJ=pfft_wave(jj,JAXIS)+partialK
        else
           partialJ=0.0
        end if
        ii=0
        n=kOffset+jOffset+1
        do i= pfft_localLimits(LOW,IAXIS),pfft_localLimits(HIGH,IAXIS)
           ii=ii+1
           partial=pfft_wave(ii,IAXIS)+partialJ
           if(partial/=0.0)partial= -1.0/partial
           
           if (((pfft_ndim==1 .AND. jOffset==0 .AND. i==0) .OR. &
                (pfft_ndim>1 .AND. jj==1 .AND. j==0))  .AND. &
               isPacked                  ) then
              
                 transArray(n) = transArray(n)  * factorArray(n)
                 transArray(n+1)=transArray(n+1)* factorArray(n+1)
              n=n+2

           else if(isComplex)then
           
              if (conjugate) then
                 transArray(n) = transArray(n)  * factorArray(n) &
                                +transArray(n+1)* factorArray(n+1)
              else
                 transArray(n) = transArray(n)  * factorArray(n) &
                                -transArray(n+1)* factorArray(n+1)
              end if
              if (conjugate) then
                 transArray(n+1)=- transArray(n)  * factorArray(n+1) &
                                  +transArray(n+1)* factorArray(n)
              else
                 transArray(n+1)=transArray(n)  * factorArray(n+1) &
                                +transArray(n+1)* factorArray(n)
              end if
              n=n+2
           else
              transArray(n)=transArray(n)*factorArray(n)
              n=n+1  
              
           end if
        end do
        jOffset=jOffset+jfactor
     end do
     kOffset=kOffset+kfactor
  end do
  
end subroutine gr_pfftMult
