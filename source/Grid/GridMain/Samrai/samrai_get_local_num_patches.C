
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
#include "Patch.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "tbox/PIO.h"
#include "PatchHierarchy.h"

// FLASH Vars
#include "Flash.h"
#include "f_to_c_names.h"

using namespace SAMRAI;



extern "C" void FTOC(samrai_get_local_num_patches)(int* num_patches)
{

  // returns number of patches on local processor
  /*
  *num_patches = 0;

  int temp_num_patches = 0;

  for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ln++) {
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getNumberOfPatches(temp_num_patches);

    *num_patches = *num_patches + temp_num_patches;

  } 
  */  

}
