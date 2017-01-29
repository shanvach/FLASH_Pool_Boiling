
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
#include "PatchHierarchy.h"
#include "Index.h"
#include "tbox/Pointer.h"
#include "tbox/PIO.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"
#include "VariableDatabase.h"


// FLASH Vars
#include "Flash.h"
#include "f_to_c_names.h"

using namespace SAMRAI;
using namespace hier;


extern "C" void FTOC(samrai_get_list_of_blocks)(int* level_no,
                      int* patch_no,
                      int* count)
{

  /*
  *count = 0;

   for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ln++) {
      Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
      for (PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
         Pointer<Patch<NDIM> > patch = level->getPatch(ip());


         int pn = patch->getPatchNumber();

       //fill up the patch and level array to be returned to f90 routine
       patch_no[*count] = pn;
       level_no[*count] = ln;

       (*count) ++;

      }
   }
  */
 
}



