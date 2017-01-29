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

double Time, Redshift, OmegaM, OmegaL, BSize, Hubble;

void FTOC(sim_readgadgetheader)(int* success,
				int* snapshot_number,
				int* files,
				char* path,
				char* basename,
				int* sim_numParticles,
				double* sim_boxSize,
				double* sim_time,
				double* sim_redshift,
				double* sim_hubble,
				double* sim_omegaMatter,
				double* sim_omegaLambda)
{

  char input_fname[200];
  int  i;

  sprintf(input_fname, "%s/%s_%03d", path, basename, (*snapshot_number));
  status = read_header(input_fname, (*files));
  
  if (status >= 0) {

    *sim_numParticles = NumPart;
    *sim_boxSize = BSize;
    *sim_time = Time;
    *sim_redshift = Redshift;
    *sim_hubble = Hubble;
    *sim_omegaMatter = OmegaM;
    *sim_omegaLambda = OmegaL;

  }

  *success = status;

  return;

}

int read_header(char *fname, int files)
{
  FILE *fd;
  char   buf[200];
  int    i,j,k,dummy,ntot_withmasses;
  int    t,n,off,pc,pc_new,pc_sph;

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

      fclose(fd);

    }

  Time = header1.time;
  Redshift = header1.time;
  
  return 0;
}











