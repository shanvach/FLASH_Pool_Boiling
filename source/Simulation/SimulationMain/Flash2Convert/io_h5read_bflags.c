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

/****if* source/IO/IOMain/hdf5/io_h5read_bflags_
 *
 * NAME
 *
 *  io_h5read_bflags_
 *
 *
 * SYNOPSIS
 *
 *  void io_h5read_bflags_(file_identifier, bflags, local_blocks, total_blocks,
 *                        global_offset)
 *
 *  void io_h5read_bflags_(hid_t *, int [][MFLAGS], int *, int *, int *)
 *
 *
 * DESCRIPTION
 *
 *  Read in a single processors's portion of the the bflags array
 *  (bflags) from a FLASH HDF5 checkpoint file.  Given the HDF5 file handle 
 *  (file_identifier), the block to start reading at (global_offset) read in
 *  local_blocks blocks worth of data.  This dataset is only necessary in
 *  Paramesh3
 *
 *
 * ARGUMENTS
 *
 *  file_identifier          the HDF5 file handle for the current file (as
 *                           returned by io_h5open_file_for_read_)
 *
 *  bflags                   arry of int flags which can be used to control computation
 *                           on the grid blocks which are inherited by children from their
 *                           parents (returned)
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

void FTOC(io_h5read_bflags)(hid_t* file_identifier,
               int bflags[][MFLAGS],
               int* local_blocks,
               int* total_blocks,
               int* global_offset)
{
  hid_t dataspace, dataset, memspace;
  herr_t status;

  int rank;
  hsize_t dimens_2d[2];

  hsize_t start_2d[2];
  hsize_t stride_2d[2], count_2d[2];

  /* open the dataset */
  dataset = H5Dopen(*file_identifier, "bflags");

  dataspace = H5Dget_space(dataset);

  /* create the hyperslab -- this will differ on the different processors */
  start_2d[0] = (hsize_t) (*global_offset);
  start_2d[1] = 0;

  stride_2d[0] = 1;
  stride_2d[1] = 1;

  count_2d[0] = (hsize_t) (*local_blocks);
  count_2d[1] = MFLAGS;

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_2d, 
                        stride_2d, count_2d, NULL);

  if (status < 0){
    printf("Error: Unable to select hyperslab for bflags\n");
    Driver_abortFlashC("Error: Unable to select hyperslab for bflags\n");
  }


  /* create the memory space */
  rank = 2;
  dimens_2d[0] = *local_blocks;
  dimens_2d[1] = MFLAGS;

  memspace = H5Screate_simple(rank, dimens_2d, NULL);

  /* read the data */
  status = H5Dread(dataset, H5T_NATIVE_INT, memspace, dataspace, 
                H5P_DEFAULT, bflags);

  if (status < 0){
    printf("Error: Unable to read dataset for bflags\n");
    Driver_abortFlashC("Error: Unable to read dataset for bflags\n");
  }

  H5Sclose(memspace); 
  H5Sclose(dataspace);
  H5Dclose(dataset);

}
