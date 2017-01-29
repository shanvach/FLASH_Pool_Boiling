#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mangle_names.h"

struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;

int NumPart, Ngas;

struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;

  float  Rho, U, Temp, Ne;
} *P;

int *Id;

void FTOC(pt_readgadgetparticles)(int* success,
				  int* snapshot_number,
				  int* files,
				  char* path,
				  char* basename,
				  double* posx,
				  double* posy,
				  double* posz,
				  double* velx,
				  double* vely,
				  double* velz,
				  double* mass,
				  double* type,
				  double* tags)
{

  char input_fname[200];
  int  i, status;

  sprintf(input_fname, "%s/%s_%03d", path, basename, (*snapshot_number));
  status = load_snapshot(input_fname, (*files));

  if (status >= 0) {

    for (i = 0; i < NumPart; i++) {

      posx[i] = (double)P[i].Pos[0];
      posy[i] = (double)P[i].Pos[1];
      posz[i] = (double)P[i].Pos[2];
    
      velx[i] = (double)P[i].Vel[0];
      vely[i] = (double)P[i].Vel[1];
      velz[i] = (double)P[i].vel[2];
      
      mass[i] = (double)P.Mass;
      
      tags[i] = (double)Id[i];

      type[i] = (double)P[i].Type;

    }

    free_memory();

  }

  *success = status;

  return;

}

/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.)
 */
int load_snapshot(char *fname, int files)
{
  FILE *fd;
  char   buf[200];
  int    i,j,k,dummy,ntot_withmasses;
  int    t,n,off,pc,pc_new,pc_sph;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i=0, pc=1; i<files; i++, pc=pc_new)
    {
      if(files>1)
	sprintf(buf,"%s.%d",fname,i);
      else
	sprintf(buf,"%s",fname);

      if(!(fd=fopen(buf,"r")))
	{
	  printf("can't open file `%s`\n",buf);
	  return -1;
	}

      printf("reading `%s' ...\n",buf); fflush(stdout);

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      if(files==1)
	{
	  for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	    NumPart+= header1.npart[k];
	  Ngas= header1.npart[0];
	}
      else
	{
	  for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	    NumPart+= header1.npartTotal[k];
	  Ngas= header1.npartTotal[0];
	}

      for(k=0, ntot_withmasses=0; k<5; k++)
	{
	  if(header1.mass[k]==0)
	    ntot_withmasses+= header1.npart[k];
	}

      if(i==0)
	allocate_memory();

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;
    

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&Id[pc_new], sizeof(int), 1, fd);
	      pc_new++;
	    }
	}
      SKIP;


      if(ntot_withmasses>0)
	SKIP;
      for(k=0, pc_new=pc; k<6; k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      P[pc_new].Type=k;

	      if(header1.mass[k]==0)
		fread(&P[pc_new].Mass, sizeof(float), 1, fd);
	      else
		P[pc_new].Mass= header1.mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses>0)
	SKIP;
      

      if(header1.npart[0]>0)
	{
	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].U, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  if(header1.flag_cooling)
	    {
	      SKIP;
	      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		{
		  fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
		P[pc_sph].Ne= 1.0;
		pc_sph++;
	      }
	}

      fclose(fd);
    }
  
  return 0;
}

/* this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(void)
{
  printf("Allocating memory from snapshot...\n");

  if(!(P=malloc(NumPart*sizeof(struct particle_data))))
    {
      fprintf(stderr,"Failed to allocate memory.\n");
      return -1;
    }
  
  P--;   /* start with offset 1 */

  
  if(!(Id=malloc(NumPart*sizeof(int))))
    {
      fprintf(stderr,"Failed to allocate memory.\n");
      return -1;
    }
  
  Id--;   /* start with offset 1 */

  printf("Allocating memory from snapshot...done\n");

  return 0;
}

/* this routine frees the memory for the 
 * particle data.
 */
int free_memory(void)
{
  printf("Freeing memory from snapshot...\n");

  free(P);
  free(Id);
  
  printf("Freeing memory from snapshot...done\n");

  return 0;
}


  











