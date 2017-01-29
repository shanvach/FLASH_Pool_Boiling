
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
!  
!  Gets the dx/dy/dz for a given blockId on the Grid
!  dx is the size of one zone in the x direction of a block
!
!  
! ARGUMENTS 
!
!  patch_no - SAMRAI patch number
!  level_no - SAMRAI level number
!  del - array of size MDIM returned holding the dx, dy, and dz values
!
*/

extern "C" void FTOC(samrai_get_deltas)( const int* patch_no,
                               const int* level_no,
                               double* del
                               )
{

  // returns deltas (dx/dy/dz) for a given patch

}
