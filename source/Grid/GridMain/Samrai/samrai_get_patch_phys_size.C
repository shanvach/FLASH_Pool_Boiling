
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

using namespace SAMRAI;

/*
! DESCRIPTION 
!  Gets the physical size of the patch along each dimension by calculating
!  the deltas of each direction times the number of zones.  It is important
!  to calculate the size of a patch consistently because rounding errors
!  can occur in some cases.  Use this function to get the size of the 
!  patch rather than calculating a new way.  Calculating the physical
!  size of the patch using deltas also is general enough to work 
!  when the patch sizes are not uniform
!  
! ARGUMENTS
!  patch_no - integer, SAMRAI patch number
!  level_no - integer, SAMRAI level number
!  patchSize(MDIM) - returned array of size MDIM holding the size of 
!              each dimension of the block
!
*/

extern "C" void FTOC(samrai_get_patch_phys_size)( const int* patch_no,
                                   const int* level_no,
                                   double* patchSize
                                   )
{

  // given a patch_no and a level_no get the physical size of patch

  //patchSize(1) = 2.0
  //patchSize(2) = 3.0
  //patchSize(3) = 4.0

 
}
