

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

!!   DESCRIPTION
!!    This subroutine is an accessor function that gets the coordinates of
!!    the cells in a given block.
!!    Coordinates are retrieved one axis at a time, meaning you can get the i, j, _or_ k
!!    coordinates with one call.  If you want all the coordinates, all axises, you
!!    need to call Grid_getCellCoords 3 times, one for each axis.
!!
!!
!! ARGUMENTS
!!            
!!   axis - specifies the integer index coordinates of the cells being retrieved.
!!          axis can have one of three different values, IAXIS, JAXIS or KAXIS 
!!          (defined in constants.h as 1,2 and 3)
!!
!!   patch_no - integer patch_number
!! 
!!   level_no - integer level number
!!
!!   edge - integer value with one of four values, LEFT, RIGHT, CENTER or WILDCARD
!!          (These are #defined constants located in constants.h)
!!        The edge argument specifies what side of the zone to get, the CENTER
!!        point, the LEFT side or the RIGHT side of the zone.  A user can get all
!!        three edges, LEFT, CENTER and RIGHT coordinates of a zone with the
!!        WILDCARD flag.
!!
!!
!!   guardcells - logical value, if true, coordinates of a block including guardcells
!!                are returned, if false, only interior cell coords are returned
!!  
!!   coordinates - The array holding the data returning the coordinate values
!!             coordinates is of size (see below)
!!           
!!   size : a 2d integer array specifying the dimensions for coordinates
!!        
!!          The values of the 'size' array depend on other arguments.
!!   
!!         size(1) = holds the number of edges that are returned either a 1 or 3
!!                 1 for 1 edge, 3 if edge = WILDCARD
!!         size(2) = size of the vector returned.  If guardcell true (and axis =IAXIS) 
!!               then
!!                   size(2) = nguard*2 + nxb for a fixed block size example.
!!                   if guardcell = false then size(2) = nxb
!***
*/

extern "C" void FTOC(samrai_get_cell_coords)( const int* axis,
                                    const int* patch_no,
                                    const int* level_no,
                                    const int* edge,
                                    const bool* guardcell,
                                    double* coordinates,
                                    const int* size
                                    )
{

  //Some of these arguments might not be needed with Samrai, basic idea is 
  // given a patch_no, level_no, and an axis, get the coordinates.


  //



}
