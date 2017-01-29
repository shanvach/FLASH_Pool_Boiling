
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
#include "VariableDatabase.h"
#include "PatchHierarchy.h"

// FLASH Vars
#include "Flash.h"
#include "f_to_c_names.h"

using namespace SAMRAI;

/*
! DESCRIPTION 
!  
!  Fills guard cells for all patches
!  
! ARGUMENTS 
! 
!  myPE         - local processor number, not sure if Samrai needs this
!  var_id     - fill GC for only this variable, 
!               if ALLVARS, fill guard cells for all variables
!               (ALLVARS is a constant defined in constants.h)
!   
!  mapperType - specifies the type of mapper, specifically if the physicaldata
!               needs conversion to another set of variables (primitive to conservative)
!
!  direction - direction of guardcell fill.  User can specify ALLDIR for all (x,y,z)
!         directions, or if for example the algorithm only does one directional
!         sweep at a time then time can be saved by filling only the guardcell
!         direction that is needed
!
!
! NOTES
!
!
!
*/

extern "C" void FTOC(samrai_fill_guard_cells)( const int* myPE,
                                     const int* var_id,
                                     const int* mapperType,
                                     const int* direction)
                                     
{
  // fills up guard cells for each patch at all levels of refinements
}
