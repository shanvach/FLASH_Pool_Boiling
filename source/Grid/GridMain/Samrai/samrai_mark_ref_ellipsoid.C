
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



extern "C" void FTOC(samrai_mark_ref_ellipsoid)( const double* ic,
                                     const double* jc,
                                     const double* kc,
                                     const double* a1,
                                     const double* a2,
                                     const double* a3,
                                     const int* lref
                                     )
{

  // Refine all patches containing points on an ellipsoidal surface centered on
  // (ic,jc,kc) with semimajor axes (a1,a2,a3).  Either patches are brought up to
  // a specific level of refinement or each patch is refined once.

}
