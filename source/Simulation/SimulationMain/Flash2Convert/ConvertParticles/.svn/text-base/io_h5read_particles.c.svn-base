
/* This file contains the functions that read the data from the HDF5 file
 * The functions accept the PARAMESH data through arguments, since C cannot
 * handle common blocks 
 */


#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "mangle_names.h"
#include <hdf5.h>
#include "hdf5_flash.h"
/*#include <mpi.h>*/
#include "Flash.h"
#include "constants.h"

typedef struct particle {
  int intProperty[2];
  double realProperty[NPART_PROPS-2];
} particle_type;

int Driver_abortFlashC(char* message);

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

void FTOC(io_get_numparticles)(hid_t *file_identifier,
		       int *numParticles)
{
  
  /* 
   * get_numParticles: look at the particles entry and get the size
   *                of the dataset based on the dimensions of the array 
   *
   * input:  file handle          (integer)
   *
   * output: number of particles  (integer)
   *
   * return value: -1 indicates an error getting the dimension
   *
   */

  herr_t   status;

  hsize_t maximum_dims[10];
  hsize_t dataspace_dims[10];

  hid_t dataset, dataspace;

  herr_t (*old_func)(void*);
  void *old_client_data;
   
  /* temporarily turn off error handling, to probe for the file's existence */
  H5Eget_auto(&old_func, &old_client_data);
  H5Eset_auto(NULL, NULL);

  /* get the dataset ID for the coordinates record, and the dataspace in
     this record */
  dataset = H5Dopen(*file_identifier, "tracer particles");

  /* restore the error handling */
  H5Eset_auto(old_func, old_client_data);

  if (dataset >= 0) {
    dataspace = H5Dget_space(dataset);

    /* get the dimensions of the dataspace */
    status = H5Sget_simple_extent_dims(dataspace, dataspace_dims, maximum_dims);
    if(status < 0){
      printf("No particles found in tracer particles dataset!\n");
      Driver_abortFlashC("No particles found in tracer particles dataset!\n");
    }

    *numParticles = dataspace_dims[0];

    H5Dclose(dataset);  
    H5Sclose(dataspace);

  } else {
    *numParticles = 0;
  }

}

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

void FTOC(io_h5read_f2_particles)(hid_t* file_identifier,
                        int* totalparticles,
                        int* localnp,
                        int* particle_offset,
                        char intPropNames[][24],
                        char realPropNames[][24],
                        particle_type particles[])
{


  hid_t    dataspace, dataset, memspace;
  herr_t   status;

  int      rank;
  hsize_t  dimens_1d, maxdims_1d;

  hssize_t start_1d;
  hsize_t  stride_1d, count_1d;

/*the particle data type*/
  hid_t    particle_tid;


/* open the dataset */

  dataset = H5Dopen(*file_identifier, "tracer particles");


/* get the particle data space */

  dataspace  = H5Dget_space(dataset);


/* get the particle data type */

  particle_tid = H5Dget_type(dataset);


/* select this processor's hyperslab */

  rank       = 1;
  maxdims_1d = (hsize_t) (*totalparticles);
  start_1d  = (hssize_t) (*particle_offset);
  stride_1d = 1;
  count_1d  = (hsize_t) (*localnp);
  dimens_1d = (hsize_t) (*localnp);

  status     = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start_1d,
                                   &stride_1d, &count_1d, NULL);

  if (status < 0){
    printf("Error: Unable to select hyperslab for tracer particles dataspace\n");
    Driver_abortFlashC("Error: Unable to select hyperslab for tracer particles dataspace\n");
  }

  memspace   = H5Screate_simple(rank, &dimens_1d, &maxdims_1d);


/* read data from the dataset */
 
  status = H5Dread(dataset, particle_tid, memspace, dataspace, H5P_DEFAULT,
                   particles);

  if (status < 0){
    printf("Error: Unable to read particles from dataset\n");
    Driver_abortFlashC("Error: Unable to read particles from dataset\n");
  }


/* release resources */

  status = H5Tclose(particle_tid);
  status = H5Sclose(memspace);
  status = H5Sclose(dataspace);
  status = H5Dclose(dataset);

}

/********************************************************************/

void FTOC(h5_read_f2_numprops)(hid_t* file_identifier, int* numintprops,
                               int* numrealprops){


  hid_t    dataspace, dataset; 

/*the particle data type*/
  hid_t    particle_tid;

  int numProps; /*overall number of properties*/
  int i;

  herr_t status;

  /*numIntProps  = (int)0;*/
  /*numRealProps = (int)0;*/

/* open the dataset */

  dataset = H5Dopen(*file_identifier, "tracer particles");


/* get the particle data space */

  dataspace  = H5Dget_space(dataset);


/* get the particle data type */

  particle_tid = H5Dget_type(dataset);

  numProps = H5Tget_nmembers(particle_tid);

  *numintprops = 0;
  *numrealprops = 0;
  
  for(i = 0; i < numProps; ++i){

    if(H5Tequal(H5Tget_member_type(particle_tid, (unsigned) i), H5T_NATIVE_INT))
      ++(*numintprops); 
    else
      ++(*numrealprops); 
  }
  
  
  if ((*numintprops)+(*numrealprops) != numProps){
    printf("Error: Unable to read particles from dataset\n");
    Driver_abortFlashC("Error: Unable to read particles from dataset\n");
  }
  
  
  /* release resources */
  
  status = H5Tclose(particle_tid);
  status = H5Sclose(dataspace);
  status = H5Dclose(dataset);
 
  return; 
}
                        
