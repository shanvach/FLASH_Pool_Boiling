#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "mangle_names.h"
#include <hdf5.h>
#include "hdf5_flash.h"
#include "Flash.h"
#include "constants.h"


int Driver_abortFlashC(char* message);

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/****if* source/IO/IOMain/hdf5/parallel/io_h5read_which_child_
 *
 * NAME
 *
 *  io_h5read_which_child_
 *
 *
 * SYNOPSIS
 *
 *  void io_h5read_which_child_(file_identifier, which_child, local_blocks, total_blocks,
 *                        global_offset)
 *
 *  void io_h5read_which_child_(hid_t *, int [], int *, int *, int *)
 *
 *
 * DESCRIPTION
 *
 *  Read in a single processors's portion of the refinement level data 
 *  (which_child) from a FLASH HDF5 checkpoint file.  Given the HDF5 file handle 
 *  (file_identifier), the block to start reading at (global_offset) read in
 *  local_blocks blocks worth of data.  
 *
 *
 * ARGUMENTS
 *
 *  file_identifier          the HDF5 file handle for the current file (as
 *                           returned by io_h5open_file_for_read_)
 *
 *  which_child              integer array identifying which part of the parents
 *                           volume nthis child corresponds to.
 *
 *  local_blocks             the number of blocks worth of data to read
 *
 *  total_blocks             the total number of blocks stored in the file
 *                           (as returned from io_h5readHeaderInfo_)
 *
 *  global_offset            the starting position of the current processor's
 *                           data (specified in blocks)
 *
 ****/

void FTOC(io_h5read_which_child)(hid_t* file_identifier,
                  int which_child[],
                  int* local_blocks,
                  int* total_blocks,
                  int* global_offset)
{
  hid_t dataspace, dataset, memspace;
  herr_t status;

  int rank;
  hsize_t dimens_1d;

  hsize_t start_1d;
  hsize_t stride_1d, count_1d;

  /* open the dataset */
  dataset = H5Dopen(*file_identifier, "which child"); 

  dataspace = H5Dget_space(dataset);

  
  /* create the hyperslab -- this will differ on the different processors */
  start_1d = (hsize_t) (*global_offset);
  stride_1d = 1;
  count_1d = (hsize_t) (*local_blocks);

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start_1d, 
                        &stride_1d, &count_1d, NULL);

  if (status < 0){
    printf("Error: Unable to select hyperslab for which_child dataspace\n");
    Driver_abortFlashC("Error: Unable to select hyperslab for which_child dataspace\n");
  }

  /* create the memory space */
  rank = 1;
  dimens_1d = *local_blocks;
  memspace = H5Screate_simple(rank, &dimens_1d, NULL);

  /* read the data */
  status = H5Dread(dataset, H5T_NATIVE_INT, memspace, dataspace, 
                H5P_DEFAULT, which_child);

  if (status < 0){
    printf("Error: Unable to read data for which_child\n");
    Driver_abortFlashC("Error: Unable to read data for which_child\n");
  }

  
  H5Sclose(memspace); 
  H5Sclose(dataspace);
  H5Dclose(dataset);

}
