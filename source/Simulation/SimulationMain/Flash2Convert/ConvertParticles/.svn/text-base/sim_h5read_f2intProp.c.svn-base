
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "mangle_names.h"
#include <hdf5.h>
#include "hdf5_flash.h"
#include "Flash.h"
#include "constants.h"

int Driver_abortFlashC(char *message);

void FTOC(sim_h5read_f2intprop)(hid_t *file_identifier,
                                char f2name[24],
                                int *localNumParts,
                                int *offset,
                                double *particles){
  
  hid_t dataspace, dataset, memspace;
  herr_t status;

  int rank,i;
  hsize_t dimens_1d, maxdimens_1d;
  hsize_t start_1d, stride_1d, count_1d;

  hid_t particle_tid;
  hid_t integer_tid;

  int *readInInts;
  char f2Cname[25];

  strncpy(f2Cname,f2name, 24);

  for(i = 0; i < 24; ++i){
    if(f2Cname[i] == ' ')
      f2Cname[i] = '\0';
  }
  printf("*********%s**********\n", f2Cname);

  readInInts = (int*)malloc(sizeof(int)*(*localNumParts));
  
  dataset = H5Dopen(*file_identifier, "tracer particles");
  particle_tid = H5Dget_type(dataset);
  dataspace = H5Dget_space(dataset);

  /*Closest thing I could find for how to get an individual member.*/

  integer_tid = H5Tcreate(H5T_COMPOUND, sizeof(int));
  status = H5Tinsert(integer_tid, f2Cname, 0, H5T_NATIVE_INT);

  start_1d = (hsize_t)(*offset);
  stride_1d = 1;
  count_1d = (hsize_t)(*localNumParts);

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start_1d,
                               &stride_1d, &count_1d, NULL);
  
  dimens_1d = (hsize_t)(*localNumParts);
  memspace = H5Screate_simple(rank, &dimens_1d, NULL);
  

  status = H5Dread(dataset, integer_tid, memspace, dataspace, H5P_DEFAULT, 
                   readInInts);

  for(i = 0; i < *localNumParts; ++i){
    particles[i] = (double)readInInts[i];
  }

  for(i = 0; i < *localNumParts; ++i){
    printf("-- %f --\n", particles[i]);
  }
  
  if(status != 0){
    printf("[IO ERROR]: Error reading in %s\n", f2Cname);
    Driver_abortFlashC("Error reading in integer particle property!\n");
  }


  free(readInInts);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Sclose(memspace);

  return;
}
