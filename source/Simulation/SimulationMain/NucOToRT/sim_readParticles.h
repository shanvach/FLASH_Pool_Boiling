#ifndef __READPARTICLES_H
#define __READPARTICLES_H
#define MASTERPE 0

/*error codes*/
#define NORMAL_STATUS 0
#define DATASET_OPEN_FAIL -1
#define FILE_OPEN_FAIL -2
#define HYPERSLAB_SELECT_FAIL -3
#define MEMSPACE_SELECT_FAIL -4
#define DATASPACE_SELECT_FAIL -5
#define DATA_READ_ERROR -6

#define PROP_STRING_LENGTH 24




//int readParticleData(hid_t file_identifier, int particle_start, 
//		     int particle_count, int npart_props, double* particles);
//int getParticleFileInfo(hid_t file_id, int *num_particles, int *num_props);
//int getPropertyNames(hid_t file_id, char *prop_names);
//int HDF5FileOpen(char* filename, hid_t* file_identifier, int create);
//void HDF5FileClose(hid_t* file_identifier);

#endif
