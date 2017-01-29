## Lines starting with ## are comments inside template file
## All other lines including empty lines are non-comments
## 
## This file is a template for generating the setup_getFlashUnits.F90
## source file.  For syntax of this file see "Readme.template".
##
##
## VALID VARIABLE NAMES FOR THIS TEMPLATE
##
## unit_names -> list of strings containing names of units 
##
## each unitname must have length atmost 80 chars 
##
!!****f* object/setup_flashUnits
!!
!! NAME
!!
!!  setup_getFlashUnits
!!
!!
!! SYNOPSIS
!!
!!
!!  setup_getFlashUnits(unit_names)
!!
!!  setup_getFlashUnits(character())
!!
!!
!! DESCRIPTION
!!
!!  Return a character array of size NUM_UNITS containing
!!  the names of all of the FLASH units used to assemble
!!  the current executable
!!
!!  The unit_names variable should be declared as
!!
!!    use flashUnits
!!
!!  
!!    character (len=MAX_STRING_LENGTH) :: flash_units(NUM_UNITS) 
!!
!!
!!  The length of each character string is set to MAX_STRING_LENGTH,
!!  which is defined in the automatically generated flash_defines.fh
!!
!!***

  subroutine setup_getFlashUnits(unit_names)

#include "constants.h"
    implicit none

    integer, PARAMETER :: NUM_UNITS = %(COUNT_unit_names)s
    character (len=MAX_STRING_LENGTH) :: unit_names(NUM_UNITS)
    character (len=80), PARAMETER :: all_unit_names(NUM_UNITS) = (/&
    &"%(unit_names !",&\n    &")s"/)

!!  if MAX_STRING_LENGTH < 80, then only the first MAX_STRING_LENGTH
!!  chars get copied over. So there is no problem with overflow here

    unit_names = all_unit_names 

    return

  end subroutine setup_getFlashUnits

  subroutine setup_getNumFlashUnits(numUnits)

#include "constants.h"
    implicit none

    integer, intent(out) :: numUnits
    integer, PARAMETER :: NUM_UNITS = %(COUNT_unit_names)s

    numUnits = NUM_UNITS

    return

  end subroutine setup_getNumFlashUnits

