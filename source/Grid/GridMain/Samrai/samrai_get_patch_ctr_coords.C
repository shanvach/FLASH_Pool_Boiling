
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

// Headers for major algorithm/data structure objects
#include "BergerRigoutsos.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CartesianVizamraiDataWriter.h"
#include "GriddingAlgorithm.h"
#include "StandardTagAndInitialize.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "LoadBalancer.h"

// FLASH Vars
#include "Flash.h"
#include "f_to_c_names.h"

using namespace SAMRAI;

/*
! DESCRIPTION 
!   Gets the coordinates of the center of the patch identified by
!   patch_no and level_no.  Returns the coordinates in an array patchCenter
!
!
! ARGUMENTS
!  patch_no - SAMRAI patch number
!  level_no - SAMRAI level number
!  patchCenter - returned array of size MDIM holding the blockCenter coords
!
*/

extern "C" void FTOC(samrai_get_patch_ctr_coords)( const int* patch_no,
                                       const int* level_no,
                                       int* patchCenter 
                                       )
{

  // returns coordinates of the center of the patch

}
