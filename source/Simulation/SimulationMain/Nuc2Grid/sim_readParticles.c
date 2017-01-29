#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "hdf5.h"
#include <string.h>
#include "sim_readParticles.h"
#include "mangle_names.h"



/*extracts the number of particles contained within a file as well as the 
  number of properties*/
void FTOC(getparticlefileinfo)(hid_t *file_id, int *num_particles, int *num_props){
  
  hid_t dataset, dataspace;
  hsize_t  dims[3];
  
  dataset = H5Dopen(*file_id, "tracer particles");
  if(dataset < 0){
    printf("Could not open tracer particles dataset!\n");
    return; 
  }
  
  dataspace = H5Dget_space(dataset);
    if(dataset < 0){
    printf("Could not open dataspace for particles !\n");
    return;
  }


  H5Sget_simple_extent_dims(dataspace, dims, NULL);

  *num_particles = (int)dims[0];
  *num_props = (int)dims[1];

  H5Sclose(dataspace);
  H5Dclose(dataset);

  return;
}

/* read in the particle property names from a file, names are 
   assumed to have a length of 24 characters.  Make sure you allocate
   enough room for the property names*/
void FTOC(getpropertynames)(hid_t *file_id, char prop_names[][24]){
 
  hid_t dataset, dataspace, memspace, string_type;
  hsize_t dims[2], maxdims[2], start[2], stride[2], count[2];
  herr_t status;
  int i, rank;

  string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, PROP_STRING_LENGTH);

  dataset = H5Dopen(*file_id, "particle names");
  dataspace = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace, dims, maxdims);

  rank = 2;
  start[0] = 0;
  start[1] = 0;

  stride[0] = 1;
  stride[1] = 1;
  
  count[0] = dims[0];
  count[1] = dims[1];
  
  memspace = H5Screate_simple(rank, dims, NULL);

  status = H5Dread(dataset, string_type, memspace,dataspace, H5P_DEFAULT, prop_names);
  if (status < 0){
    Driver_abortFlashC("Could not read in particle names from nuc file!\n");
  }
  // for(i=0; i<dims[0]; ++i){
    // strncat(prop_names[i], "                       ", 24);
  // printf("prop_names: %s\n", prop_names[i]);
  //}
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  H5Tclose(string_type);


  return;
}


/*opens up a file.  Takes a filename, a file identifier, and a 
  flag indicating whether or not we are creating a given file.
  1 = create, 0 = just open. The file is opened in parallel across
  all processors with this.*/

void FTOC(hdf5fileopen)(char filename[81], hid_t* file_identifier, int *create){

    int ierr;
    MPI_Info fileInfoTemplate;
    hid_t acc_template;

    /*set up file access*/
    ierr = MPI_Info_create(&fileInfoTemplate);

    acc_template = H5Pcreate(H5P_FILE_ACCESS);
    /*open for a parallel read*/

    ierr = H5Pset_fapl_mpio(acc_template, MPI_COMM_WORLD, fileInfoTemplate);

      if(!(*create))
	*file_identifier = H5Fopen(strtok(filename, " "), H5F_ACC_RDWR, acc_template);
    else
      *file_identifier = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, acc_template);
    
    ierr = H5Pclose(acc_template);
    
    if(*file_identifier < 0){
        fprintf(stderr, "Error opening file, %s!\nError returned is %d!\n",
	         filename, *file_identifier);
	Driver_abortFlashC("could not open file %s\n",filename);
	return;
    }
    
    return;
}


/*close a hdf5 file out.*/
void FTOC(hdf5fileclose)(hid_t* file_identifier){
  
  
  H5Fclose(*file_identifier);
  
  return;
}

