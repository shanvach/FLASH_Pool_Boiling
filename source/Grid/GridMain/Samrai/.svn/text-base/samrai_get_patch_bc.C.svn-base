
#include "SAMRAI_config.h"

#include <stdlib.h>
#include <string>
#include <fstream>
using namespace std;

#include <sys/stat.h>

// Headers for basic SAMRAI objects
#include "tbox/SAMRAIManager.h"
#include "tbox/Database.h"
#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "tbox/MPI.h"
#include "Patch.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "tbox/PIO.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"
#include "VariableDatabase.h"



// FLASH Vars
#include "Flash.h"
#include "f_to_c_names.h"
#include "constants.h"

using namespace SAMRAI;

/*
! DESCRIPTION 
!  Returns the boundary condition for each face of the patch.
!  
!  this function finds out if a patch face is
!  on the physical boundary, and if so, returns the boundary condition.
!  The boundary conditions are defined in the header file
!  constants.h, ie. OUTFLOW, REFLECTING, PERIODIC.  NOT_BOUNDARY, also
!  defined in constants.h is returned if a patch face is not on a 
!  physical boundary.
!   
! ARGUMENTS 
!
!  patch_no - SAMRAI patch Number
!  level_no - SAMRAI level Number
!  faces(2,MDIM)   - array returned holding boundary conditions
!                
!            the first index of the array can take on values LOWERFACE or
!            UPPERFACE, and the second index can be IAXIS, JAXIS or KAXIS
!
! NOTES
!
!   The #define constants LOWERFACE, UPPERFACE, IAXIS, JAXIS and KAXIS
!   are defined in constants.h and are
!   meant to ease the readability of the code.  
!   instead of faces(2,3) = PERIODIC the code reads
!   faces(UPPERFACE, KAXIS) = PERIODIC
!
*/
extern "C" void FTOC(samrai_get_patch_bc)( const int* patch_no,
                                 const int* level_no,
                                 int* faces)
                                 
{

  // given a patch_no and level_no, figures out whether the patch 
  // is on the physical boundary 
  // and if so returns boundary condition 

  //example
  //faces(LOWERFACE, IAXIS) = PERIODIC
  //faces(UPPERFACE, IAXIS) = NOT_BOUNDARY 
  //faces(LOWERFACE, JAXIS) = PERIODIC
  //faces(UPPERFACE, JAXIS) = NOT_BOUNDARY 
  //faces(LOWERFACE, KAXIS) = NOT_BOUNDARY
  //faces(UPPERFACE, KAXIS) = NOT_BOUNDARY 
  
}
