
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
!  
!  Gets simulation data in storage for specified cells and variables
!  
! ARGUMENTS 
!
!  axis(MDIM) !! specifies the integer index coordinates of the cells
!                A WILDCARD along any of the axes tell the funtion
!                to get all the data along that dimension. 
!                If the axis doesn't exist (for example KAXIS in a 2D 
!                problem, it value defaults to KHI_GC, which should be 1.
!  
!  patch_no : SAMRAI patch number
!
!  level_no : SAMRAI level number
!
!  guarcell : is true if guardcell data is to be included, otherwise false
!
!  variable : Specifies which variable to get. A WILDCARD gets all variables.
!
!  datablock : real array containing the data, 
!             the size defined by the parameter "size"
!
!  size : (maybe call something else) the array dimension specs for datablock
!          the valid sizes for all dimensions are :
!          IF (variable == WILDCARD) size(1)=NVAR ELSE size(1)=1 
!          IF (axis(1)==WILDCARD || guardcell==true) size(2)=IHI_GC 
!          IF (axis(1)==WILDCARD || guardcell==false) size(2)=NXB 
!          IF (axis(1)/=WILDCARD) size(2) = 1
!          IF (axis(2)==WILDCARD || guardcell==true) size(3)=JHI_GC 
!          IF (axis(2)==WILDCARD || guardcell==false) size(3)=NYB 
!          IF (axis(2)/=WILDCARD) size(3) = 1
!          IF (axis(3)==WILDCARD || guardcell==true) size(4)=KHI_GC 
!          IF (axis(3)==WILDCARD || guardcell==false) size(4)=NZB 
!          IF (axis(3)/=WILDCARD) size(4) = 1
!
!        Significance:  
!         IF only one axis is WILDCARD, fetches one row of data 
!         IF any two axes are WILDCARD, fetches one Plane of data
!         IF all three axes = WILDCARD, fetches a complete Block
! 
! EXAMPLES :For example <WILDCARD,2,KHI_GC> in axis
!           for a 2D problem tells the subroutine to fetch data of cells (row)
!           in the second column from the interior of the block if guardcell
!           is false, otherwise it fetches in the second column of guardcells. 
!           Similarly <WILDCARD,WILDCARD,5> in 3D fetches the 5th XY plane data.
!           from the interior if guardcell is false, or 5-NGUARD-th plane
!           if guardcell is true.
!           If all the values in axis are WILDCARDS, 
!           then the entire block data is returned.
! 
!
*/

extern "C" void FTOC(samrai_get_data)( const int* axis,
                               const int* patch_no,
                               const int* level_no,
                               const bool* guardcell,
                               const int* variable,
                               double* dataBlock,
                               const int* size
                               )
{

  // returns data in storage, for specified patchs' cells and variables

}
