#include "hdf5.h"
#include "mpi.h"
#include <ctype.h>

#include "underscore.h"
#define float double
#define H5T_NATIVE_FLOAT H5T_NATIVE_DOUBLE

void CREATE_File_matlab_r8 (char* file_name_in, hid_t* file_id) {

  hid_t plist_id;
  char file_name[80];
  int i; int ii;

  /*
     Don't know why I had to add this code, but it appears that some
     non-numeric, non-chararcter data gets inserted into this file name
     when called from the c_interface.  This code will remove any weird
     characters
  */
  for(i=0; i<80; i++) {
    if (!isalnum(file_name_in[i]) && !ispunct(file_name_in[i])) {
      for (ii=i;ii<79;ii++) {
        file_name_in[ii]=file_name_in[ii+1];
      }
    }
  }

  /* strip off fortran junk from the filename */
  strncpy(file_name, file_name_in, 80);
  for(i=79; i>=0; i--) {
    if (file_name[i] != ' ') { file_name[i+1]=0; break; }
  }

  /* Set up file access property list with parllel I/O access */
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

  /* Create a new file collectively */
  *file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

}



void WRITE_Header_matlab_r8(hid_t* file_id, int* tot_blocks,
			    int* ndim,
			    int* nxb, int* nyb, int* nzb,
                            int* nvar,
                            int* nvar_fc,
                            int* nvar_ec,
                            int* nvar_nc,
			    float* user_attr_1,
			    float* user_attr_2,
			    float* user_attr_3,
			    float* user_attr_4,
			    float* user_attr_5) {

  hsize_t  dims_header[1];
  hid_t    dset_header, dataspace_header, group_id, attribute_id;
  herr_t   status;
  int      mpi_rank;
  int      temp[9];

  dims_header[0] = 9;
  dataspace_header = H5Screate_simple(1, dims_header, NULL);
  dset_header = H5Dcreate(*file_id, "Global Data", H5T_NATIVE_INT, 
			  dataspace_header, H5P_DEFAULT);

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  if (mpi_rank == 0) {
    temp[0] = *tot_blocks;
    temp[1] = *ndim;
    temp[2] = *nxb;
    temp[3] = *nyb;
    temp[4] = *nzb;
    temp[5] = *nvar;
    temp[6] = *nvar_fc;
    temp[7] = *nvar_ec;
    temp[8] = *nvar_nc;
    status = H5Dwrite(dset_header, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
                      H5P_DEFAULT, temp);
  }

  H5Sclose(dataspace_header);
  H5Dclose(dset_header);

  /* ADD USER ATTRIBUTES */
  group_id = H5Gopen(*file_id, "/");

  /* attribute 1 */
  dataspace_header = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "user attribute 1", H5T_NATIVE_FLOAT,
			   dataspace_header, H5P_DEFAULT);
  if (mpi_rank == 0) {
    status = H5Awrite(attribute_id, H5T_NATIVE_FLOAT, user_attr_1);
  }
  H5Sclose(dataspace_header);
  H5Aclose(attribute_id);

  /* attribute 2 */
  dataspace_header = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "user attribute 2", H5T_NATIVE_FLOAT,
			   dataspace_header, H5P_DEFAULT);
  if (mpi_rank == 0) {
    status = H5Awrite(attribute_id, H5T_NATIVE_FLOAT, user_attr_2);
  }
  H5Sclose(dataspace_header);
  H5Aclose(attribute_id);

  /* attribute 3 */
  dataspace_header = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "user attribute 3", H5T_NATIVE_FLOAT,
			   dataspace_header, H5P_DEFAULT);
  if (mpi_rank == 0) {
    status = H5Awrite(attribute_id, H5T_NATIVE_FLOAT, user_attr_3);
  }
  H5Sclose(dataspace_header);
  H5Aclose(attribute_id);

  /* attribute 4 */
  dataspace_header = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "user attribute 4", H5T_NATIVE_FLOAT,
			   dataspace_header, H5P_DEFAULT);
  if (mpi_rank == 0) {
    status = H5Awrite(attribute_id, H5T_NATIVE_FLOAT, user_attr_4);
  }
  H5Sclose(dataspace_header);
  H5Aclose(attribute_id);

  /* attribute 5 */
  dataspace_header = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "user attribute 5", H5T_NATIVE_FLOAT,
			   dataspace_header, H5P_DEFAULT);
  if (mpi_rank == 0) {
    status = H5Awrite(attribute_id, H5T_NATIVE_FLOAT, user_attr_5);
  }
  H5Sclose(dataspace_header);
  H5Aclose(attribute_id);

  H5Gclose(group_id);
      
}

 
void CREATE_Block_data_type_matlab_r8(hid_t* block_type_id, 
				   int* mdim, 
				   int* ndim, 
				   int* ngid, 
				   int* mflags, 
				   int* nvar_chk_cc, 
				   int* nvar_chk_fc, 
				   int* nvar_chk_ec, 
				   int* nvar_chk_nc, 
				   int* il0, int* iu0, 
				   int* jl0, int* ju0, 
				   int* kl0, int* ku0) 
{
  
  int rank_array;
  hsize_t dims_array[1];
  int k2d, k3d;

  hid_t gid_tid, bflags_tid, coord_tid, bnd_box_tid, unk_tid;
  hid_t facevarx_tid, facevary_tid, facevarz_tid;
  hid_t unk_e_x_tid, unk_e_y_tid, unk_e_z_tid;
  hid_t unk_n_tid;
  int nx, ny, nz;
  int size_block_type;
  int offset;
  
  
  k2d = 0;
  if (*ndim >= 2) k2d = 1;
  k3d = 0;
  if (*ndim == 3) k3d = 1;

  /* Create the memory data type */

  size_block_type = (sizeof(float)*3) + 
    (sizeof(float) * *ngid) * + (sizeof(float) * *mflags) +
    (sizeof(float) * *mdim) + (sizeof(float) * 2 * *mdim) +
    (sizeof(float)) + 
    (sizeof(float) * *nvar_chk_cc * (*iu0-*il0) * (*ju0-*jl0) * (*ku0-*kl0)) +
    (sizeof(float) * *nvar_chk_fc * (*iu0-*il0+1) * (*ju0-*jl0) * (*ku0-*kl0)) +
    (sizeof(float) * *nvar_chk_fc * (*iu0-*il0) * (*ju0-*jl0+k2d) * (*ku0-*kl0)) +
    (sizeof(float) * *nvar_chk_fc * (*iu0-*il0) * (*ju0-*jl0) * (*ku0-*kl0+k3d)) +
    (sizeof(float) * *nvar_chk_ec * (*iu0-*il0) * (*ju0-*jl0+k2d) * (*ku0-*kl0+k3d)) +
    (sizeof(float) * *nvar_chk_ec * (*iu0-*il0+1) * (*ju0-*jl0) * (*ku0-*kl0+k3d)) +
    (sizeof(float) * *nvar_chk_ec * (*iu0-*il0+1) * (*ju0-*jl0+k2d) * (*ku0-*kl0)) +
    (sizeof(float) * *nvar_chk_nc * (*iu0-*il0+1) * (*ju0-*jl0+k2d) * (*ku0-*kl0+k3d));

  *block_type_id = H5Tcreate (H5T_COMPOUND, size_block_type);

  rank_array = 1;
  dims_array[0] = *ngid;
  gid_tid = H5Tarray_create(H5T_NATIVE_FLOAT, rank_array, 
			    dims_array, NULL);
  rank_array = 1;
  dims_array[0] = *mflags;
  bflags_tid = H5Tarray_create(H5T_NATIVE_FLOAT, rank_array, 
			       dims_array, NULL);
  rank_array = 1;
  dims_array[0] = *mdim;
  coord_tid = H5Tarray_create(H5T_NATIVE_FLOAT, rank_array, 
                              dims_array, NULL);
  rank_array = 1;
  dims_array[0] = 2 * *mdim;
  bnd_box_tid = H5Tarray_create(H5T_NATIVE_FLOAT, rank_array, 
				dims_array, NULL);
  rank_array = 1;
  nx = *iu0 - *il0;
  ny = *ju0 - *jl0;
  nz = *ku0 - *kl0;

  if (*nvar_chk_cc > 0) {
  dims_array[0] = *nvar_chk_cc * nx * ny * nz;
  unk_tid = H5Tarray_create(H5T_NATIVE_FLOAT, rank_array, 
			    dims_array, NULL);
  }

  if (*nvar_chk_fc > 0) {
  dims_array[0] = *nvar_chk_fc * (nx+1) * ny * nz;
  facevarx_tid = H5Tarray_create(H5T_NATIVE_FLOAT, rank_array, 
				 dims_array, NULL);
  dims_array[0] = *nvar_chk_fc * nx * (ny+k2d) * nz;
  facevary_tid = H5Tarray_create(H5T_NATIVE_FLOAT, rank_array, 
				 dims_array, NULL);
  dims_array[0] = *nvar_chk_fc * nx * ny * (nz+k3d);
  facevarz_tid = H5Tarray_create(H5T_NATIVE_FLOAT, rank_array, 
				 dims_array, NULL);
  }

  if (*nvar_chk_ec > 0) {
  dims_array[0] = *nvar_chk_ec * nx * (ny+k2d) * (nz+k3d);
  unk_e_x_tid = H5Tarray_create(H5T_NATIVE_FLOAT, rank_array, 
				 dims_array, NULL);
  dims_array[0] = *nvar_chk_ec * (nx+1) * ny * (nz+k3d);
  unk_e_y_tid = H5Tarray_create(H5T_NATIVE_FLOAT, rank_array, 
				 dims_array, NULL);
  dims_array[0] = *nvar_chk_ec * (nx+1) * (ny+k2d) * nz;
  unk_e_z_tid = H5Tarray_create(H5T_NATIVE_FLOAT, rank_array, 
				 dims_array, NULL);
  }

  if (*nvar_chk_nc > 0) {
  dims_array[0] = *nvar_chk_nc * (nx+1) * (ny+k2d) * (nz+k3d);
  unk_n_tid = H5Tarray_create(H5T_NATIVE_FLOAT, rank_array, 
				 dims_array, NULL);
  }

  offset = 0;

  H5Tinsert(*block_type_id, "lrefine",  
	    offset, H5T_NATIVE_FLOAT);
  offset = offset + sizeof(float);

  H5Tinsert(*block_type_id, "nodetype", 
	    offset, H5T_NATIVE_FLOAT);
  offset = offset + sizeof(float);

  H5Tinsert(*block_type_id, "which_child", 
	    offset, H5T_NATIVE_FLOAT);
  offset = offset + sizeof(float);

  H5Tinsert(*block_type_id, "gid", 
	    offset, gid_tid);
  offset = offset + (sizeof(float) * *ngid);

  H5Tinsert(*block_type_id, "bflags",  
	    offset, bflags_tid);
  offset = offset + (sizeof(float) * *mflags);

  H5Tinsert(*block_type_id, "coord", 
	    offset, coord_tid);
  offset = offset + (sizeof(float) * *mdim);

  H5Tinsert(*block_type_id, "bnd_box", 
	    offset , bnd_box_tid);
  offset = offset + (sizeof(float) * 2 * *mdim);

  H5Tinsert(*block_type_id, "work_block", 
	    offset, H5T_NATIVE_FLOAT);
  offset = offset + sizeof(float);

  if (*nvar_chk_cc > 0) {
    H5Tinsert(*block_type_id, "unk", 
	      offset, unk_tid);
    offset = offset + (sizeof(float) * *nvar_chk_cc * (*iu0-*il0) * (*ju0-*jl0) * (*ku0-*kl0));
  }
  if (*nvar_chk_fc > 0) {
    H5Tinsert(*block_type_id, "facevarx", 
	      offset, facevarx_tid);
    offset = offset + (sizeof(float) * *nvar_chk_fc * (*iu0-*il0+1) * (*ju0-*jl0) * (*ku0-*kl0));
    H5Tinsert(*block_type_id, "facevary", 
	      offset, facevary_tid);
    offset = offset + (sizeof(float) * *nvar_chk_fc * (*iu0-*il0) * (*ju0-*jl0+k2d) * (*ku0-*kl0));
    H5Tinsert(*block_type_id, "facevarz", 
	      offset, facevarz_tid);
    offset = offset + (sizeof(float) * *nvar_chk_fc * (*iu0-*il0) * (*ju0-*jl0) * (*ku0-*kl0+k3d));
  }
  if (*nvar_chk_ec > 0) {
    H5Tinsert(*block_type_id, "unk_e_x", 
	      offset, unk_e_x_tid);
    offset = offset + (sizeof(float) * *nvar_chk_ec * (*iu0-*il0) * (*ju0-*jl0+k2d) * (*ku0-*kl0+k3d));
    H5Tinsert(*block_type_id, "unk_e_y", 
	      offset, unk_e_y_tid);
    offset = offset + (sizeof(float) * *nvar_chk_ec * (*iu0-*il0+1) * (*ju0-*jl0) * (*ku0-*kl0+k3d));
    H5Tinsert(*block_type_id, "unk_e_z", 
	      offset, unk_e_z_tid);
    offset = offset + (sizeof(float) * *nvar_chk_ec * (*iu0-*il0+1) * (*ju0-*jl0+k2d) * (*ku0-*kl0));
  }
  if (*nvar_chk_nc > 0) {
    H5Tinsert(*block_type_id, "unk_n", 
	      offset, unk_n_tid);
  }

  H5Tclose(gid_tid);
  H5Tclose(bflags_tid);
  H5Tclose(coord_tid);
  H5Tclose(bnd_box_tid);
  if (*nvar_chk_cc > 0) {
    H5Tclose(unk_tid);
  }
  if (*nvar_chk_fc > 0) {
    H5Tclose(facevarx_tid);
    H5Tclose(facevary_tid);
    H5Tclose(facevarz_tid);
  }
  if (*nvar_chk_ec > 0) {
    H5Tclose(unk_e_x_tid);
    H5Tclose(unk_e_y_tid);
    H5Tclose(unk_e_z_tid);
  }
  if (*nvar_chk_nc > 0) {
    H5Tclose(unk_n_tid);
  }

}


void CREATE_Block_data_set_matlab_r8(hid_t* dset_block, 
				  hid_t* file_id,
				  hid_t* block_type_id,
				  int*   tot_blocks_wr)
{

  hsize_t dims_block[1];
  hid_t   dataspace_block;

  /* Create the dataspace for block info */
  dims_block[0] = *tot_blocks_wr;
  dataspace_block = H5Screate_simple(1, dims_block, NULL);
  
  /* Create the dataset for block and cell centered data */
  /* Create a dataset collectively for each processor */
  *dset_block = H5Dcreate(*file_id, "Block Data", *block_type_id, 
                          dataspace_block, H5P_DEFAULT);
  H5Sclose(dataspace_block);

}


void WRITE_Blocks_matlab_r8(hid_t* file_id, 
			 hid_t* dset_block, 
			 hid_t* block_type_id, 
			 int* max_lnblocks, 
			 int* n_to_left, 
			 int* lnblocks, 
			 int* mdim,
			 int* ndim,
			 int* ngid,
			 int* mflags,
			 int lrefine[], 
			 int nodetype[], 
			 int which_child[], 
			 int gid[], 
			 int bflags[],
			 float coord[],
			 float bnd_box[],
			 float work_block[],
			 float unk[],
			 int*  nvar,
			 int*  nvar_chk_cc,
			 int   checkp_on_cc[],
			 float facevarx[], float facevary[], float facevarz[],
			 int*  nbndvar,
			 int*  nvar_chk_fc,
			 int   checkp_on_fc[],
			 float unk_e_x[], float unk_e_y[], float unk_e_z[],
			 int*  nbndvare,
			 int*  nvar_chk_ec,
			 int   checkp_on_ec[],
			 float unk_n[],
			 int*  nbndvarc,
			 int*  nvar_chk_nc,
			 int   checkp_on_nc[],
			 int*  nx,  int* ny, int* nz,
			 int*  il0, int* iu0,
			 int*  jl0, int* ju0,
			 int*  kl0, int* ku0)
{
  hid_t  memspace, filespace;
  hsize_t  count[1];
  hssize_t offset[1];
  hid_t plist_id;
  int ivar, i, j, k, lb, i2, j2, k2, nx2, ny2, ivar2;
  int istart, iend, jstart, jend, kstart, kend;
  herr_t status;
  int k2d, k3d;
  int off, size_block_type;

  float * block_pointer;

  k2d = 0;
  if (*ndim >= 2) k2d = 1;
  k3d = 0;
  if (*ndim == 3) k3d = 1;

  size_block_type = (sizeof(float)*3) + 
    (sizeof(float) * *ngid) * + (sizeof(float) * *mflags) +
    (sizeof(float) * *mdim) + (sizeof(float) * 2 * *mdim) +
    (sizeof(float)) + 
    (sizeof(float) * *nvar_chk_cc * (*iu0-*il0) * (*ju0-*jl0) * (*ku0-*kl0)) +
    (sizeof(float) * *nvar_chk_fc * (*iu0-*il0+1) * (*ju0-*jl0) * (*ku0-*kl0)) +
    (sizeof(float) * *nvar_chk_fc * (*iu0-*il0) * (*ju0-*jl0+k2d) * (*ku0-*kl0)) +
    (sizeof(float) * *nvar_chk_fc * (*iu0-*il0) * (*ju0-*jl0) * (*ku0-*kl0+k3d)) +
    (sizeof(float) * *nvar_chk_ec * (*iu0-*il0) * (*ju0-*jl0+k2d) * (*ku0-*kl0+k3d)) +
    (sizeof(float) * *nvar_chk_ec * (*iu0-*il0+1) * (*ju0-*jl0) * (*ku0-*kl0+k3d)) +
    (sizeof(float) * *nvar_chk_ec * (*iu0-*il0+1) * (*ju0-*jl0+k2d) * (*ku0-*kl0)) +
    (sizeof(float) * *nvar_chk_nc * (*iu0-*il0+1) * (*ju0-*jl0+k2d) * (*ku0-*kl0+k3d));

  block_pointer = (float *) malloc(size_block_type);
  
  /* Create property list for collective dataset write */
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  /*H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);*/
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

  nx2 = *iu0 - *il0;
  ny2 = *ju0 - *jl0;
  kstart = *kl0;
  kend   = *ku0;
  jstart = *jl0;
  jend   = *ju0;
  istart = *il0;
  iend   = *iu0;
  if (*ndim < 3) {
    kstart = 0;
    kend   = 1;
  }
  if (*ndim < 2) {
    kstart = 0;
    kend   = 1;
    jstart = 0;
    jend   = 1;
  }
  
  for (lb = 0; lb < *max_lnblocks; lb++) {

    if (lb < *lnblocks) {
      off = 0;
      *block_pointer = (float)lrefine[lb];
      off = off + 1;
      *(block_pointer+off) = (float)nodetype[lb];
      off = off + 1;
      *(block_pointer+off) = (float)which_child[lb];
      off = off + 1; 
      for (i = 0; i < *ngid; i++) {
	*(block_pointer+off) = (float)gid[(*ngid * lb) + i];
	off = off + 1; 
      }
      for (i = 0; i < *mflags; i++) {
	*(block_pointer+off) = (float)bflags[(*mflags * lb) + i];
	off = off + 1;
      }
      for (i = 0; i < *mdim; i++) {
	*(block_pointer+off) = coord[(*mdim * lb) + i];
	off = off + 1;
      }
      for (i = 0; i < 2 * *mdim; i++) {
	*(block_pointer+off) = bnd_box[(2 * *mdim * lb) + i];
	off = off + 1;
      }
      *(block_pointer+off) = work_block[lb];
      off = off + 1;

      /* unk data */

      if (*nvar_chk_cc > 0) {

      k2 = 0;
      for (k = kstart; k < kend; k++) {
	j2 = 0;
	for (j = jstart; j < jend; j++) {
	  i2 = 0;
	  for (i = istart; i < iend; i++) {
	    ivar2 = 0;
	    for (ivar=0; ivar < *nvar; ivar++) {

	      if (checkp_on_cc[ivar] == 1) {
		*(block_pointer+off) =
		  unk[(*nz * *ny * *nx * *nvar * lb) + 
		      (*ny * *nx * *nvar * k) + 
		      (*nx * *nvar * j) + 
		      (*nvar * i) + 
		      ivar];
		ivar2 = ivar2 + 1;
		off = off + 1;
	      }

	    } /* end ivar loop */
	    i2 = i2 + 1;
	  } /* end i loop */
	  j2 = j2 + 1;
	} /* end j loop */
	k2 = k2 + 1;
      } /* end k loop */ 

      } /* end if (nvar_chk_cc */

      /* facevarx data */

      if (*nvar_chk_fc > 0) {

      k2 = 0;
      for (k = kstart; k < kend; k++) {
	j2 = 0;
	for (j = jstart; j < jend; j++) {
	  i2 = 0;
	  for (i = istart; i < iend+1; i++) {
	    ivar2 = 0;
	    for (ivar = 0; ivar < *nbndvar; ivar++) {

	      if (checkp_on_fc[ivar + 2 * ivar] == 1) {
                *(block_pointer+off) =
		  facevarx[(*nz * *ny * (*nx+1) * *nbndvar * lb) +
			   (*ny * (*nx+1) * *nbndvar * k) +
			   ((*nx+1) * *nbndvar * j) +
			   (*nbndvar * i) +
			   ivar];
		ivar2 = ivar2 + 1;
		off = off + 1;
	      }

	    } /* end ivar loop */
	    i2 = i2 + 1;
	  } /* end i loop */
	  j2 = j2 + 1;
	} /* end j loop */
	k2 = k2 + 1;
      } /* end k loop */

      /* facevary data */

      k2 = 0;
      for (k = kstart; k < kend; k++) {
	j2 = 0;
	for (j = jstart; j < jend+k2d; j++) {
	  i2 = 0;
	  for (i = istart; i < iend; i++) {
	    ivar2 = 0;
	    for (ivar = 0; ivar < *nbndvar; ivar++) {

	      if (checkp_on_fc[ivar + 1 + 2 * ivar] == 1) {
                *(block_pointer+off) =
		  facevary[(*nz * (*ny+k2d) * *nx * *nbndvar * lb) +
			   ((*ny+k2d) * *nx * *nbndvar * k) +
			   (*nx * *nbndvar * j) +
			   (*nbndvar * i) +
			   ivar];
		ivar2 = ivar2 + 1;
		off = off + 1;
	      }
	    } /* end ivar loop */
	    i2 = i2 + 1;
	  } /* end i loop */
	  j2 = j2 + 1;
	} /* end j loop */
	k2 = k2 + 1;
      } /* end k loop */

      /* facevarz data */

      k2 = 0;
      for (k = kstart; k < kend+k3d; k++) {
	j2 = 0;
	for (j = jstart; j < jend; j++) {
	  i2 = 0;
	  for (i = istart; i < iend; i++) {
	    ivar2 = 0;
	    for (ivar = 0; ivar < *nbndvar; ivar++) {

	      if (checkp_on_fc[ivar + 2 + 2 * ivar] == 1) {
                *(block_pointer+off) =
		  facevarz[((*nz+k3d) * *ny * *nx * *nbndvar * lb) +
			   (*ny * *nx * *nbndvar * k) +
			   (*nx * *nbndvar * j) +
			   (*nbndvar * i) +
			   ivar];
		ivar2 = ivar2 + 1;
		off = off + 1;
	      }

	    } /* end ivar loop */
	    i2 = i2 + 1;
	  } /* end i loop */
	  j2 = j2 + 1;
	} /* end j loop */
	k2 = k2 + 1;
      } /* end k loop */

      } /* end if (nvar_chk_fc */

      /* unk_e_x data */

      if (*nvar_chk_ec > 0) {

      k2 = 0;
      for (k = kstart; k < kend+k3d; k++) {
	j2 = 0;
	for (j = jstart; j < jend+k2d; j++) {
	  i2 = 0;
	  for (i = istart; i < iend; i++) {
	    ivar2 = 0;
	    for (ivar = 0; ivar < *nbndvare; ivar++) {

	      if (checkp_on_ec[ivar + 2 * ivar] == 1) {
                *(block_pointer+off) =
		  unk_e_x[((*nz+k3d) * (*ny+k2d) * *nx * *nbndvare * lb) +
			  ((*ny+k2d) * *nx * *nbndvare * k) +
			  (*nx * *nbndvare * j) +
			  (*nbndvare * i) +
			  ivar];
		ivar2 = ivar2 + 1;
		off = off + 1;
	      }

	    } /* end ivar loop */
	    i2 = i2 + 1;
	  } /* end i loop */
	  j2 = j2 + 1;
	} /* end j loop */
	k2 = k2 + 1;
      } /* end k loop */

      /* unk_e_y data */

      k2 = 0;
      for (k = kstart; k < kend+k3d; k++) {
	j2 = 0;
	for (j = jstart; j < jend; j++) {
	  i2 = 0;
	  for (i = istart; i < iend+1; i++) {
	    ivar2 = 0;
	    for (ivar = 0; ivar < *nbndvare; ivar++) {

	      if (checkp_on_ec[ivar + 1 + 2 * ivar] == 1) {
                *(block_pointer+off) =
		  unk_e_y[((*nz+k3d) * *ny * (*nx+1) * *nbndvare * lb) +
			  (*ny * (*nx+1) * *nbndvare * k) +
			  ((*nx+1) * *nbndvare * j) +
			  (*nbndvare * i) +
			  ivar];
		ivar2 = ivar2 + 1;
		off = off + 1;
	      }

	    } /* end ivar loop */
	    i2 = i2 + 1;
	  } /* end i loop */
	  j2 = j2 + 1;
	} /* end j loop */
	k2 = k2 + 1;
      } /* end k loop */

      /* unk_e_z data */

      k2 = 0;
      for (k = kstart; k < kend; k++) {
	j2 = 0;
	for (j = jstart; j < jend+k2d; j++) {
	  i2 = 0;
	  for (i = istart; i < iend+1; i++) {
	    ivar2 = 0;
	    for (ivar = 0; ivar < *nbndvare; ivar++) {

	      if (checkp_on_ec[ivar + 2 + 2 * ivar] == 1) {
                *(block_pointer+off) =
		  unk_e_z[(*nz * (*ny+k2d) * (*nx+1) * *nbndvare * lb) +
			  ((*ny+k2d) * (*nx+1) * *nbndvare * k) +
			  ((*nx+1) * *nbndvare * j) +
			  (*nbndvare * i) +
			  ivar];
		ivar2 = ivar2 + 1;
		off = off + 1;
	      }

	    } /* end ivar loop */
	    i2 = i2 + 1;
	  } /* end i loop */
	  j2 = j2 + 1;
	} /* end j loop */
	k2 = k2 + 1;
      } /* end k loop */

      } /* end if (nvar_chk_ec */

      /* unk_n data */

      if (*nvar_chk_nc > 0) {

      k2 = 0;
      for (k = kstart; k < kend+k3d; k++) {
	j2 = 0;
	for (j = jstart; j < jend+k2d; j++) {
	  i2 = 0;
	  for (i = istart; i < iend+1; i++) {
	    ivar2 = 0;
	    for (ivar=0; ivar < *nbndvarc; ivar++) {

	      if (checkp_on_nc[ivar] == 1) {
		*(block_pointer+off) =
		  unk_n[((*nz+k3d) * (*ny+k2d) * (*nx+1) * *nbndvarc * lb) + 
			((*ny+k2d) * (*nx+1) * *nbndvarc * k) + 
			((*nx+1) * *nbndvarc * j) + 
			(*nbndvarc * i) + 
			ivar];
		ivar2 = ivar2 + 1;
		off = off + 1;
	      }

	    } /* end ivar loop */
	    i2 = i2 + 1;
	  } /* end i loop */
	  j2 = j2 + 1;
	} /* end j loop */
	k2 = k2 + 1;
      } /* end k loop */

      } /* end if (nvar_chk_nc */

      count[0] = 1;
      memspace = H5Screate_simple(1, count, NULL);
  
      offset[0] = *n_to_left + lb; 
      filespace = H5Dget_space(*dset_block);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, 
                          NULL, count, NULL);

    } else {

      int lb2;

      lb2 = *lnblocks - 1;

      count[0] = 1;
      memspace = H5Screate_simple(1, count, NULL);
  
      offset[0] = *n_to_left + lb2; 
      filespace = H5Dget_space(*dset_block);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, 
                          NULL, count, NULL);

    }
    
    status = H5Dwrite(*dset_block, *block_type_id, memspace, filespace, 
                      plist_id, block_pointer);

    H5Sclose(filespace);
    H5Sclose(memspace);
    
  }

  H5Pclose(plist_id);

  free(block_pointer);

}

#ifdef UNDERSCORE
void write_blocks_matlab_r8_ (char* file_name_in, 
#endif
#ifdef DOUBLE_UNDERSCORE
void write_blocks_matlab_r8__ (char* file_name_in, 
#endif
#ifndef UNDERSCORE 
#ifndef DOUBLE_UNDERSCORE
void write_blocks_matlab_r8 (char* file_name_in, 
#endif
#endif
			   int*  tot_blocks,
			   int*  tot_blocks_wr,
			   int*  max_lnblocks,
			   int*  lnblocks,
			   int*  n_to_left,
			   int*  mdim,
			   int*  ndim,
			   int*  ngid,
			   int*  mflags,
			   int   lrefine[], 
			   int   nodetype[], 
			   int   which_child[], 
			   int   gid[], 
			   int   bflags[],
			   float coord[],
			   float bnd_box[],
			   float work_block[],
			   float unk[],
			   int*  nvar,
			   int*  nvar_chk_cc,
			   int   checkp_on_cc[],
			   float facevarx[], float facevary[], float facevarz[],
			   int*  nbndvar,
			   int*  nvar_chk_fc,
			   int   checkp_on_fc[],
			   float unk_e_x[], float unk_e_y[], float unk_e_z[],
			   int*  nbndvare,
			   int*  nvar_chk_ec,
			   int   checkp_on_ec[],
			   float unk_n[],
			   int*  nbndvarc,
			   int*  nvar_chk_nc,
			   int   checkp_on_nc[],
			   int*  nx,  int* ny, int* nz,
			   int*  il0, int* iu0,
			   int*  jl0, int* ju0, 
			   int*  kl0, int* ku0,
			   float* user_attr_1,
			   float* user_attr_2,
			   float* user_attr_3,
			   float* user_attr_4,
			   float* user_attr_5)
{
  hid_t  file_id, dset_block;  /* file and dataset ids */
  hid_t  block_type_id;

  CREATE_File_matlab_r8(file_name_in, &file_id);

  WRITE_Header_matlab_r8(&file_id, tot_blocks, ndim,
                         nx, ny, nz, 
			 nvar_chk_cc,
                         nvar_chk_fc,
                         nvar_chk_ec,
                         nvar_chk_nc,
			 user_attr_1,
			 user_attr_2,
			 user_attr_3,
			 user_attr_4,
			 user_attr_5
			 );

  CREATE_Block_data_type_matlab_r8(&block_type_id, 
				 mdim, 
				 ndim, 
				 ngid, 
				 mflags,
				 nvar_chk_cc, 
				 nvar_chk_fc, 
				 nvar_chk_ec, 
				 nvar_chk_nc, 
				 il0, iu0, 
				 jl0, ju0, 
				 kl0, ku0); 

  CREATE_Block_data_set_matlab_r8(&dset_block, 
				&file_id, 
				&block_type_id, 
				tot_blocks_wr);
  
  WRITE_Blocks_matlab_r8(&file_id, 
		       &dset_block, 
		       &block_type_id, 
		       max_lnblocks, 
		       n_to_left, 
		       lnblocks, 
		       mdim, 
		       ndim, 
		       ngid, 
		       mflags, 
		       lrefine, 
		       nodetype, 
		       which_child, 
		       gid, 
		       bflags, 
		       coord, 
		       bnd_box, 
		       work_block, 
		       unk,
		       nvar, 
		       nvar_chk_cc,
		       checkp_on_cc,
		       facevarx, facevary, facevarz,
		       nbndvar, 
		       nvar_chk_fc,
		       checkp_on_fc,
		       unk_e_x, unk_e_y, unk_e_z,
		       nbndvare, 
		       nvar_chk_ec,
		       checkp_on_ec,
		       unk_n,
		       nbndvarc, 
		       nvar_chk_nc,
		       checkp_on_nc,
		       nx, ny, nz,
		       il0, iu0, jl0, ju0, kl0, ku0);
  
  /* Close resources */
  H5Dclose(dset_block);
  H5Fclose(file_id);
  H5Tclose(block_type_id);

}



     
