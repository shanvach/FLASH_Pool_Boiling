
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

void FTOC(sim_h5read_f2partpropnames)(hid_t *file_identifier,
                                     char intPropNames[][24],
                                     char realPropNames[][24]){




  hid_t    dataspace, dataset, memspace;
  herr_t   status;

  int      rank;
  hsize_t  dimens_1d, maxdims_1d;

  hssize_t start_1d;
  hsize_t  stride_1d, count_1d;

  /*the particle data type*/
  hid_t    particle_tid;

  int intPropCount, realPropCount;
  int numMembers;
  char *name;
  int i,j;
  

  name = NULL;
 
  dataset = H5Dopen(*file_identifier, "tracer particles");
  dataspace = H5Dget_space(dataset);
  /*get data type*/
  particle_tid = H5Dget_type(dataset);

  numMembers = H5Tget_nmembers(particle_tid);
  
  intPropCount = realPropCount = 0;

  for(i=0; i < numMembers; i++){
      
    name = H5Tget_member_name(particle_tid, (unsigned int)i);
  
    if(H5Tequal(H5Tget_member_type(particle_tid, i), H5T_NATIVE_INT))
      strncpy(intPropNames[intPropCount++],name,24);
    else if(H5Tequal(H5Tget_member_type(particle_tid, i), H5T_NATIVE_DOUBLE))
      strncpy(realPropNames[realPropCount++], name, 24);
    else
      Driver_abortFlashC("Error: Types containted within the particle data structure are not native to this architecture.  Cannot convert particle data at this time.");
    

    free(name);
  
  }
  /*Fortranify the names.  Instead of nulls, we need whitespace*/
  
  for (i=0; i<intPropCount; ++i){
    for(j=0; j<24; ++j){
      if(intPropNames[i][j] == '\0')
        intPropNames[i][j] = ' ';
    }
  }
   
  
  for (i=0; i<realPropCount; ++i){
    for(j=0; j<24; ++j){
      if(realPropNames[i][j] == '\0')
        realPropNames[i][j] = ' ';
    }
  }
  
  H5Sclose(dataspace);
  H5Dclose(dataset);
  return;
}
