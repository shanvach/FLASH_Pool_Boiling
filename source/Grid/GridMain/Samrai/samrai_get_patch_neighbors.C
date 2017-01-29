
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
!  Gets information about patch's neighbors
!  
! ARGUMENTS 
!
!  patch_no - SAMRAI patch number
!  level_no - SAMRAI level_no
!
!  NB: This description may not be very accurate
!  FOR PARAMESH
!  neigh - neigh(1,:) contains the local blockId of the block
!        neigh(2,:) contains the processor it is on
!        There is provision for two neighbors along each face. If there is
!        only one neighbor along a face, then all entries for the second 
!        neighbor along the face contain a NULL value
*/


extern "C" void FTOC(samrai_get_patch_neighbors)( const int* patch_no,
                                      const int* level_no,
                                      int* neigh
                                      )
{

  
  // returns information on neighbors of the patch
  // 

}
