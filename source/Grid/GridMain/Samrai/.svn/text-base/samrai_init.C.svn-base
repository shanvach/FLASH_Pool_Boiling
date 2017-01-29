
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
#include "tbox/Array.h"
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
using namespace hier;
using namespace tbox;
using namespace mesh;
using namespace geom;

extern "C" void FTOC(samrai_init)(const char var_ids[][5],
                          const int* num_vars,
                          const int* nghosts,
                          const double* xminmax,
                          const int* iminmax,
                          const int* processor_layout,
                          const int* max_levels,
                          const int* refine_ratio,
                          const int* largest_patch_size,
                          const int* smallest_patch_size,
                          const double* efficiency_tolerance,
                          const double* combine_efficiency,
                          int* argc,
                          char** argv) {






  printf("inside of samrai_init\n");



   /*
    * 1) initialize MPI and SAMRAI manager
  
   tbox::MPI::init(&argc, char*argv[]);
   tbox::SAMRAIManager::startup();
*/

   /*
    * 2) create variables
    *
    * The argument "num_vars" is supplied.  Here, we create an array of 
    * Cell centered SAMRAI variables that store those that are used in FLASH.
    * We create the variables, assign a single "default" context, and register
    * them with SAMRAI's variable database.
    */

  /*
  Array<CellVariable<NDIM,double> > flash_vars;
  flash_vars.resizeArray(num_vars);

  char var_id_new[*num_vars][5];
  int i = 0;

 
  for (int i = 0; i < flash_vars.getSize(); i++) {
  
    //translate from fortran to C strings
    //the variable names are 4 characters long -- copy this into 
    //var_id_new, the 5th character is for the \0 termination 
    strncpy(var_id_new[i], var_ids[i],4);
    *(var_id_new[i] + 4) = '\0';

     flash_vars[i] = 
        new CellVariable<NDIM,double>(var_id_new[i],1);
  }


  // a single context for all variables
  Pointer<VariableContext> default_cxt = 
    VariableDatabase<DIM>::getDatabase()->getContext("DEFAULT");
  const IntVector<DIM> ghosts(nghosts);

  Array<int> var_ids[flash_vars.getSize()];
  for (int i = 0; i < var_ids.getSize(); i++) {
  
    var_id[i] = 
      variable_db->registerVariableAndContext(flash_vars[i],
                                    default_cxt);
  }
  */
  /*
   * 3) Create hierarchy
   *
   *   3a) Create a CartesianGridGeometry using xmin/xmax inputs.
   *   3b) Create a patch hierarchy
   *   3c) Create an error detector (where you tag cells for refinement)
   *   3d) Create a box generator - standard Berger-Rigoutsos
   *   3e) Create a load balancer
   *   3f) Create a Gridding Algorithm
   *   3g) Use gridding algorithm to build hierarchy
   */

  // 3a - create a Cartesian grid geometry
  /*Pointer<Database> cart_geom_db;

  // set domain_boxes
  Index<NDIM> domain_ilo;
  Index<NDIM> domain_ihi;
  for (i = 0; i < NDIM; i++) {
     domain_ilo(i) = iminmax[2*i];
     domain_ihi(i) = iminmax[2*i+1];
  }
  Box<NDIM> domain_box(domain_ilo,domain_ihi);
  cart_geom_db->putBox("domain_boxes",domain_box);

  // set x_lo, x_hi
  Array<double> domain_xlo(NDIM);
  Array<double> domain_xhi(NDIM);
  for (i = 0; i < NDIM; i++) {
     domain_xlo[i] = xminmax[2*i];
     domain_xhi[i] = xminmax[2*i+1];
  }
  cart_geom_db->putDoubleArray("x_lo",domain_xlo);
  cart_geom_db->putDoubleArray("x_hi",domain_xhi);

  Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry =
    new geom::CartesianGridGeometry<NDIM>("CartesianGeometry",
                                cart_geom_db);               

   // 3b - create a patch hierarchy
   Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy = 
      new hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

   // 3c - create an error detector
   Pointer<Database> stai_db;
   stai_db->addString("tagging_method","GRADIENT_DETECTOR");

   tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > error_detector =
      new mesh::StandardTagAndInitialize<NDIM>(
             "StandardTagAndInitialize",
             NULL,
             stai_db);

    // 3d - create a box generator
    tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator = 
       new mesh::BergerRigoutsos<NDIM>();

    // 3e - create a load balancer 
    Pointer<Database> lb_db;
    // if you want to be able to specify processor layout this is where you'd 
    // do it.
    // Array<int> processor_layout(NDIM);
    // <set processor_layout>
    // lb_db->putIntegerArray("processor_layout",processor_layout);
    tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer = 
       new mesh::LoadBalancer<NDIM>("LoadBalancer", 
                                    lb_db);

    // 3f - create a gridding algorithm
    Pointer<Database> gridding_db;
    gridding_db->putInteger("max_levels",max_levels);
    Pointer<Database> ratio_to_coarser_db;
    ratio_to_coarser_db->putIntegerArray("level_1",refine_ratio,NDIM);
    gridding_db->putDatabase(ratio_to_coarser_db);
    Pointer<Database> largest_patch_size_db;
    largest_patch_size_db->putIntegerArray("level_0",largest_patch_size,NDIM);
    gridding_db->putDatabase(largest_patch_size_db);
    Pointer<Database> smallest_patch_size_db;
    smallest_patch_size_db->putIntegerArray("level_0",smallest_patch_size,NDIM);
    gridding_db->putDatabase(smallest_patch_size_db);
    gridding_db->putDouble("efficiency_tolerance", efficiency_tolerance);
    gridding_db->putDouble("combine_efficiency", combine_efficiency);

    Pointer< mesh::GriddingAlgorithm<NDIM> > gridding_algorithm = 
       new mesh::GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                         gridding_db,
                                         error_detector,
                                         box_generator,
                                         load_balancer);

    // 3g - use gridding algorithm to build hierarchy
    Array<int> *tag_buffer_array = 
       new tbox::Array<int>(gridding_algorithm->getMaxLevels());
    int tag_buffer = 1; // later, we can set this in FLASH and pass it
    // as a runtime param
    for (int il = 0; il < gridding_algorithm->getMaxLevels(); il++) {
       (*tag_buffer_array)[il] = tag_buffer;
    }

    gridding_algorithm->makeCoarsestLevel(patch_hierarchy,loop_time);

    bool done = false;
    bool initial_time = true;
    for (int ln = 0; 
         gridding_algorithm->levelCanBeRefined(ln) && !done; 
         ln++) {
       gridding_algorithm->makeFinerLevel(patch_hierarchy,
                                          0.0,// time = 0.
                                          true, // initial_time = true
                                          (*tag_buffer_array)[ln]);
       done = !(patch_hierarchy->finerLevelExists(ln));
    }
  */

  printf("leaving samrai_init\n");

}
     
  




