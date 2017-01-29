
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
!  Gets the physical domain bounding box of the patch identified 
!  by patch_no and level_no.
!  The boundbox is defined as the lower left and upper right
!  corners of the patch
!
! ARGUMENTS
!  patch_no - integer, SAMRAI patch number
!  level_no - integer, SAMRAI level number
!  boundBox(2, MDIM) - returned array holding the boundBox coordinates in
!             each dimension
!
*/

extern "C" void FTOC(samrai_get_patch_bound_box)( const int* patch_no,
                                      const int* level_no,
                                      double* boundBox
                                      )
{

  //Samrai seems to do this a little differently, seems like first I need
  // to get the box object? - then get the coordinates of the boundBox
  // from that object


  //get the patch
  /*
  tbox::Pointer< hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(*level_no);
  tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(*patch_no);
  */
  //getBox gets the indicies -- we want physical physical bound box
  //not sure how to get physical coords, check the documentation



  
  //boundBox(1, IAXIS) = 0.0
  //boundBox(2, IAXIS) = 0.0
  //boundBox(1, JAXIS) = 1.0
  //boundBox(2, JAXIS) = 1.0


  //kda: don't forget indicies flip going from C to fortran

}
