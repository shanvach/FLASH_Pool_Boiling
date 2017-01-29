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
#include <mpi.h>
#include "Flash.h"
#include "constants.h"
/*#include "driverAPI.h"*/
/*
typedef struct particle {
  int intProperty[NUMINTPROPS];
  double realProperty[NUMREALPROPS];
} particle_type;
*/
int Driver_abortFlashC(char* message);
typedef struct sim_params_t{
    int total_blocks;
    int nsteps;
    int nxb;
    int nyb;
    int nzb;
    double time;
    double timestep;
    double redshift;
} sim_params_t;

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/****if* source/Simlation/SimulationMain/Flash2Convert/h5_read_header_info_
 *
 * NAME
 *
 *  io_h5read_header_
 *
 *
 * SYNOPSIS
 *
 *  void io_h5read_f2header_(MyPE, file_identifier, total_blocks, time, 
 *                            timestep, redshift, nsteps)
 *
 *  void io h5read_f2header_(int *, hid_t *, int *, double *, double *, 
 *                            double *, int *)
 *
 *
 * DESCRIPTION
 *
 *  This has been brought in from Flash 2 to read in Flash 2 header information
 *  which is very different from Flash 3 header information in hdf5 files.
 *  
 *
 *  Given an HDF5 file identifier (*file_identifier), read the FLASH header
 *  information from the file.  This header information is packed in a
 *  HDF5 compound object, which is constructed by this routine before reading.
 *
 *
 * ARGUMENTS
 *
 *  MyPE              the local processor number
 *
 *  file_identifer    the HDF5 file handle (as returned from 
 *                    h5_open_file_for_read_)
 *
 *  total_blocks      the total number of blocks found in the file (returned)
 *
 *  time              the simulation time (returned)
 *
 *  timestep          the current simulation timestep (returned)
 *
 *  redshift          the simulation redshift (returned)
 *
 *  nsteps            the number of steps taken to this point in the simulation
 *                    (returned)
 *
 *
 ****/

void FTOC(io_h5read_f2header)(int* MyPE,
			  hid_t* file_identifier,    /* file handle */
			  int* total_blocks,         /* total # of blocks */
			  double* time,              /* simulation time */
			  double* timestep,          /* current timestep */
			  double* redshift,          /* simulation redshift */
			  int* nsteps)               /* total # of timestep */            

{
  hid_t dataset;
  herr_t status;

  int file_version;

  hid_t sp_type;

  sim_params_t sim_params;

  /* file version number */
  dataset = H5Dopen(*file_identifier, "file format version");
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
		   H5P_DEFAULT, &file_version);
  if(status < 0){
    printf("Error: Unable to read file_version\n");
    Driver_abortFlashC("Error: Unable to read file_version\n");
  }

  H5Dclose(dataset);
  if (*MyPE == 0) {
    printf ("determined file format version to be %d.\n", file_version);
  }

  dataset = H5Dopen(*file_identifier, "simulation parameters");

  /* create the HDF 5 compound data type to describe the record */
  if (file_version > 2) {
    sp_type = H5Tcreate(H5T_COMPOUND, sizeof(sim_params_t));
  }
  else {
    sp_type = H5Tcreate(H5T_COMPOUND, sizeof(sim_params_t)-sizeof(double));
  }

  H5Tinsert(sp_type, 
            "total blocks", 
            offsetof(sim_params_t, total_blocks),
            H5T_NATIVE_INT);
 
  H5Tinsert(sp_type,
            "time",
            offsetof(sim_params_t, time),
            H5T_NATIVE_DOUBLE);

  H5Tinsert(sp_type,
            "timestep",
            offsetof(sim_params_t, timestep),
            H5T_NATIVE_DOUBLE);

  if (file_version > 2) {
    H5Tinsert(sp_type,
              "redshift",
              offsetof(sim_params_t, redshift),
              H5T_NATIVE_DOUBLE);
  }
  
  H5Tinsert(sp_type,
            "number of steps",
            offsetof(sim_params_t, nsteps),
            H5T_NATIVE_INT);

  H5Tinsert(sp_type,
            "nxb",
            offsetof(sim_params_t, nxb),
            H5T_NATIVE_INT);

  H5Tinsert(sp_type,
            "nyb",
            offsetof(sim_params_t, nyb),
            H5T_NATIVE_INT);
  
  H5Tinsert(sp_type,
            "nzb",
            offsetof(sim_params_t, nzb),
            H5T_NATIVE_INT);

  
  status = H5Dread(dataset, sp_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
		   &sim_params);

  if (status < 0) {
    printf("Error: Unable to read simulation parameters\n");
    Driver_abortFlashC("Error: Unable to read simulation parameters\n");
   
  }

  H5Tclose(sp_type);
  H5Dclose(dataset);

  *total_blocks = sim_params.total_blocks;
  *time         = sim_params.time;
  *timestep     = sim_params.timestep;
  *nsteps       = sim_params.nsteps;
  if (file_version > 2) {
    *redshift = sim_params.redshift;
  }
  else {
    *redshift = 0.0;
  }

}

