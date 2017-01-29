
#include "SAMRAI_config.h"

#include <stdlib.h>
#include <string>
#include <fstream>
using namespace std;

#include <sys/stat.h>

// Headers for basic SAMRAI objects
#include "tbox/SAMRAIManager.h"
#include "Patch.h"
#include "PatchLevel.h"
#include "PatchHierarchy.h"
#include "tbox/Pointer.h"
#include "tbox/PIO.h"
#include "VariableDatabase.h"




// FLASH Vars
#include "Flash.h"
#include "f_to_c_names.h"

// Header for application-specific algorithm/data structure object
/*#include "ConvDiff.h"*/
using namespace SAMRAI;

extern "C" void FTOC(samrai_get_patch_extent)(const int* patch_no,
                           const int* level_no,
                           int* range)  // int of size 6)
{


  // Access patch from hierarchy
  /*
  tbox::Pointer< hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(*level_no);
  tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(*patch_no);

  
  const hier::Index<NDIM> ifirst = patch->getBox().lower();
  const hier::Index<NDIM> ilast  = patch->getBox().upper();

  range[0] = ifirst(0); 
  range[1] = ilast(0); 
  range[2] = ifirst(1); 
  range[3] = ilast(1); 
  range[4] = ifirst(2); 
  range[5] = ilast(2); 
  */
 }


