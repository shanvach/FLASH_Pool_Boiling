
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
!  Apply the user-defined refinment critera to determine which patch needs 
!  to be refined and derefined.  Once the patches are marked, call 
!  amr_refine_derefine to actually carry out the refinements.  During this
!  stage, the patchess are redistributed across processors (if needed).  
!
!  After the refinement, the newly created child patches are filled via
!  prolongation from the coarse parents.  This prolongation step uses 
!  a user-defined prolongation routine (these can be picked in the Modules
!  file).
!
!  Once the prolongation is done, the guardcells are filled.  Finally, the
!  EOS is called on the patch interiors to make them thermodynamically
!  consistent.
!
!
! PARAMETERS 
!
!  nrefs - The number of steps between refinements
!  nstep - number of timesteps
!  patchList - List of all patches
!
*/


extern "C" void FTOC(samrai_update_refinement)( const int* nrefs,
                                    const int* nsteps,
                                    int* patchList
                                    )
{

  // Applies user defined criteria to determine which patches 
  // needs to be refined or derefined

}
