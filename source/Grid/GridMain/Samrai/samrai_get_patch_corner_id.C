
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
! !!DEV: verify this description.  what about guardcells 
! Returns the global integer indices of the start of the interior zone 
! of the block
! and the stride of indices along each dimensions.
! For example, in a 1 dimensional UG case with 2 blocks and nxb=8
! The index for block 1 = 1 and the index for block 2 = 9 
! 
! 
! ARGUMENTS 
!
!  patch_no - SAMRAI patch number
!  level_no - SAMRAI level number 
!  Index    - global integer indices of start of the interior zone
!             of the block
!      
!  stride   - stride for indices, in UG stride is always = 1
!           for the paramesh stride may be more than 1 depending
!             on how far down you are in the tree, In SAMRAI need to 
!             figure it out
!
*/


extern "C" void FTOC(samrai_get_patch_corner_id)( const int* patch_no,
                                      const int* level_no,
                                      int* index,
                                      int* stride)
                                      
{

  // may not be needed in samrai.  Not sure how this is implemented
  // returns global integer 'index' of the beginning of the interior of patch
  // and also the stride of indices along each dimensions

}
