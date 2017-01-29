
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

/****if* source/IO/IOMain/hdf5/io_h5read_nodetype_
 *
 * NAME
 *
 *  io_h5read_nodetype_
 *
 *
 * SYNOPSIS
 *
 *  void io_h5read_nodetype_(file_identifier, nodetype, local_blocks, total_blocks,
 *                        global_offset)
 *
 *  void io_h5read_nodetype_(hid_t *, int [], int *, int *, int *)
 *
 *
 * DESCRIPTION
 *
 *  Read in a single processors's portion of the tree node type data 
 *  (nodetype) from a FLASH HDF5 checkpoint file.  Given the HDF5 file handle 
 *  (file_identifier), the block to start reading at (global_offset) read in
 *  local_blocks blocks worth of data.  
 *
 *
 * ARGUMENTS
 *
 *  file_identifier          the HDF5 file handle for the current file (as
 *                           returned by io_h5open_file_for_read_)
 *
 *  nodetype                 the node type data (returned)
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

void FTOC(io_h5read_nodetype)(hid_t* file_identifier,
                   int nodetype[],
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
  dataset = H5Dopen(*file_identifier, "node type"); 

  dataspace = H5Dget_space(dataset);

  
  /* create the hyperslab -- this will differ on the different processors */
  start_1d = (hsize_t) (*global_offset);
  stride_1d = 1;
  count_1d = (hsize_t) (*local_blocks);

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start_1d, 
                        &stride_1d, &count_1d, NULL);

  if (status < 0){
    printf("Error: Unable to select hyperslab for nodetype\n");
    Driver_abortFlashC("Error: Unable to select hyperslab for nodetype\n");
  }


  /* create the memory space */
  rank = 1;
  dimens_1d = *local_blocks;
  memspace = H5Screate_simple(rank, &dimens_1d, NULL);

  /* read the data */
  status = H5Dread(dataset, H5T_NATIVE_INT, memspace, dataspace, 
                H5P_DEFAULT, nodetype);

  if (status < 0){
    printf("Error: Unable to read dataset for nodetype\n");
    Driver_abortFlashC("Error: Unable to read dataset for nodetype\n");
  }


  H5Sclose(memspace); 
  H5Sclose(dataspace);
  H5Dclose(dataset);

}

