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


int Driver_abortFlashC(char* message);


/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

void FTOC(io_h5read_bndbox)(hid_t* file_identifier,
                     int* maximum_blocks,   
                     double* bnd_box,
                     int* local_blocks,
                     int* total_blocks,
                     int* global_offset)
{
  hid_t dataspace, dataset, memspace;
  herr_t status;

  int rank;
  hsize_t dimens_3d[3], maxdimens_3d[3];

  hsize_t start_3d[3];
  hsize_t stride_3d[3], count_3d[3];

  int ierr;

  /* open the dataset */  
  dataset = H5Dopen(*file_identifier, "bounding box");
  
  dataspace = H5Dget_space(dataset);

  H5Sget_simple_extent_dims(dataspace, dimens_3d, maxdimens_3d);
  
  /* create the hyperslab -- this will differ on the different processors */
  start_3d[0] = (hsize_t) (*global_offset);
  start_3d[1] = 0;
  start_3d[2] = 0;

  stride_3d[0] = 1;
  stride_3d[1] = 1;
  stride_3d[2] = 1;

  count_3d[0] = (hsize_t) (*local_blocks);
  count_3d[1] = dimens_3d[1];
  count_3d[2] = dimens_3d[2];

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_3d, 
                        stride_3d, count_3d, NULL);

  if (status < 0){
    printf("Error: Unable to select hyperslab for bndbox dataspace\n");
    Driver_abortFlashC("Error: Unable to select hyperslab for bndbox dataspace\n");
  }

  /* create the memory space */
  rank = 3;
  dimens_3d[0] = *maximum_blocks;
    dimens_3d[1] = MDIM;
    dimens_3d[2] = 2;
 
  
  memspace = H5Screate_simple(rank, dimens_3d, NULL);

  start_3d[0] = 0;
  start_3d[1] = 0;
  start_3d[2] = 0;

  stride_3d[0] = 1;
  stride_3d[1] = 1;
  stride_3d[2] = 1;

  /*count_3d[0] = *local_blocks;
  /*count_3d[1] = dimens_3d[1];*/
  /*printf ("output dimensions %d\n", dimens_3d[1]);*/
  /*exit(1);
  count_3d[2] = dimens_3d[2];*/

  ierr = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
                             start_3d, stride_3d, count_3d, NULL);

  if (ierr < 0){
    printf("Error: Unable to select hyperslab for bndbox memspace\n");
    Driver_abortFlashC("Error: Unable to select hyperslab for bndbox memspace\n");
  }

  /* read the data */
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
               H5P_DEFAULT, bnd_box);

  if (status < 0){
    printf("Error: Unable to read bnd_box from dataset\n");
    Driver_abortFlashC("Error: Unable to read bnd_box from dataset\n");
  }

  H5Sclose(memspace); 
  H5Sclose(dataspace);
  H5Dclose(dataset);

}

