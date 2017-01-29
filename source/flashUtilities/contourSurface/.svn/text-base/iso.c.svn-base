/*
 * This program implements the marching cubes surface tiler described by
 * Lorensen & Cline in the Siggraph 87 Conference Proceedings.
 *
 * This program was adaptated from isovis code by Mike Krogh, NCSA, 1990.
 * Timur Linde, ASCI Flash Center, University of Chicago, 2004.
 *
*/

#include <stdio.h>
#include <math.h>
#include "cell_table.h"

#include "mangle_names.h"

/**************************** Temporary Globals **********************/

static double DATA1,DATA2,DATA3,DATA4,DATA5,DATA6,DATA7,DATA8;
static int XDIMYDIM;


/**************************** iso_surface ****************************/

void FTOC(iso_surface)(data,xdim,ydim,zdim,xmin,ymin,zmin,xmax,ymax,zmax,threshold,area)
double *data;
int *xdim,*ydim,*zdim;
double *xmin,*ymin,*zmin;
double *xmax,*ymax,*zmax;
double *threshold;
double *area;
{

  register int x,y,z;
  register int xdim1,ydim1,zdim1;
  int index;
  int npolys;
  int old_verts;
  int edge[13];
  double crossings[13][3];
  void get_cell_verts();
  void get_cell_polys();
  void calc_index_and_temps();

  double dx,dy,dz;

  dx = ((*xmax)-(*xmin))/((*xdim)-1);
  dy = ((*ymax)-(*ymin))/((*ydim)-1);
  dz = ((*zmax)-(*zmin))/((*zdim)-1);

  zdim1 = (*zdim)-1;
  ydim1 = (*ydim)-1;
  xdim1 = (*xdim)-1;

  XDIMYDIM = (*xdim)*(*ydim);

  npolys=0;
  (*area)=0.0;

  for (z=0;z<zdim1;z++)
    for (y=0;y<ydim1;y++)
      for (x=0;x<xdim1;x++) {
        calc_index_and_temps(data,x,y,z,*xdim,*ydim,*zdim,*threshold,&index);
        if (index) {
          get_cell_verts(index,x,y,z,0.,0.,0.,*threshold,crossings);
          get_cell_polys(index,&npolys,area,crossings,dx,dy,dz);
        }
      }
}


/**************************** calc_index_and_temps ********************/

void calc_index_and_temps(data,x1,y1,z1,xdim,ydim,zdim,threshold,index)
register double *data;
int x1,y1,z1;
int xdim,ydim,zdim;
register double threshold;
int *index;
{

  register double *tmp;

  *index = 0;

  tmp = data + (z1*XDIMYDIM) + (y1*xdim) + x1;

  *index += (threshold <= (DATA1 = *(tmp)));
  *index += (threshold <= (DATA2 = *(tmp + 1))) * 2;

  tmp += xdim;
  *index += (threshold <= (DATA3 = *(tmp + 1))) * 4;
  *index += (threshold <= (DATA4 = *(tmp))) * 8;

  tmp = tmp - xdim + XDIMYDIM;
  *index += (threshold <= (DATA5 = *(tmp))) * 16;
  *index += (threshold <= (DATA6 = *(tmp + 1))) * 32;

  tmp += xdim;
  *index += (threshold <= (DATA7 = *(tmp + 1))) * 64;
  *index += (threshold <= (DATA8 = *(tmp))) * 128;
 
}


/**************************** get_cell_verts ****************************/

void get_cell_verts(index,x1,y1,z1,xtrans,ytrans,ztrans,threshold,crossings)
int index;
int x1,y1,z1;
double xtrans,ytrans,ztrans;
double threshold;
double crossings[13][3];
{

  register int i;
  register int x2,y2,z2;
  int nedges;
  int crnt_edge;
 
#define linterp(a1,a2,a,b1,b2) ((double)(((a-a1) * (double)(b2-b1) / (a2-a1)) + (double)b1))

  x2 = x1+1;
  y2 = y1+1;
  z2 = z1+1;

  nedges = cell_table[index].nedges;
  for (i=0;i<nedges;i++) {
     crnt_edge = cell_table[index].edges[i];
     switch (crnt_edge) {
	case 1:
		crossings[1][0] = linterp(DATA1,DATA2,threshold,x1,x2)+xtrans;
		crossings[1][1] = (double)y1+ytrans;
		crossings[1][2] = (double)z1+ztrans;
		break;

	case 2:
		crossings[2][1] = linterp(DATA2,DATA3,threshold,y1,y2)+ytrans;
		crossings[2][0] = (double)x2+xtrans;
		crossings[2][2] = (double)z1+ztrans;
		break;

	case 3:
		crossings[3][0] = linterp(DATA4,DATA3,threshold,x1,x2)+xtrans;
		crossings[3][1] = (double)y2+ytrans;
		crossings[3][2] = (double)z1+ztrans;
		break;

	case 4:
		crossings[4][1] = linterp(DATA1,DATA4,threshold,y1,y2)+ytrans;
		crossings[4][0] = (double)x1+xtrans;
		crossings[4][2] = (double)z1+ztrans;
		break;

	case 5:
		crossings[5][0] = linterp(DATA5,DATA6,threshold,x1,x2)+xtrans;
		crossings[5][1] = (double)y1+ytrans;
		crossings[5][2] = (double)z2+ztrans;
		break;

	case 6:
		crossings[6][1] = linterp(DATA6,DATA7,threshold,y1,y2)+ytrans;
		crossings[6][0] = (double)x2+xtrans;
		crossings[6][2] = (double)z2+ztrans;
		break;

	case 7:
		crossings[7][0] = linterp(DATA8,DATA7,threshold,x1,x2)+xtrans;
		crossings[7][1] = (double)y2+ytrans;
		crossings[7][2] = (double)z2+ztrans;
		break;

	case 8:
		crossings[8][1] = linterp(DATA5,DATA8,threshold,y1,y2)+ytrans;
		crossings[8][0] = (double)x1+xtrans;
		crossings[8][2] = (double)z2+ztrans;
		break;

	case 9:
		crossings[9][2] = linterp(DATA1,DATA5,threshold,z1,z2)+ztrans;
		crossings[9][1] = (double)y1+ytrans;
		crossings[9][0] = (double)x1+xtrans;
		break;

	case 10:
		crossings[10][2] = linterp(DATA2,DATA6,threshold,z1,z2)+ztrans;
		crossings[10][1] = (double)y1+ytrans;
		crossings[10][0] = (double)x2+xtrans;
		break;

	case 11:
		crossings[11][2] = linterp(DATA4,DATA8,threshold,z1,z2)+ztrans;
		crossings[11][1] = (double)y2+ytrans;
		crossings[11][0] = (double)x1+xtrans;
		break;

	case 12:
		crossings[12][2] = linterp(DATA3,DATA7,threshold,z1,z2)+ztrans;
		crossings[12][1] = (double)y2+ytrans;
		crossings[12][0] = (double)x2+xtrans;
		break;

     }
  }
}


/**************************** get_cell_polys ****************************/

void get_cell_polys(index,npolys,area,crossings,dx,dy,dz)
int index;
int *npolys;
double *area;
double crossings[13][3];
double dx,dy,dz;
{

  register int num_o_polys;
  register int poly;
  double *p1,*p2,*p3;

  double V1x,V1y,V1z;
  double V2x,V2y,V2z;
  double Ax,Ay,Az;

  num_o_polys = cell_table[index].npolys;

  for (poly=0;poly<num_o_polys;poly++) {

    p1 = &crossings[cell_table[index].polys[(poly*3)]][0];
    p2 = &crossings[cell_table[index].polys[(poly*3)+1]][0];
    p3 = &crossings[cell_table[index].polys[(poly*3)+2]][0];

    V1x = (p2[0]-p1[0])*dx;
    V1y = (p2[1]-p1[1])*dy;
    V1z = (p2[2]-p1[2])*dz;

    V2x = (p3[0]-p1[0])*dx;
    V2y = (p3[1]-p1[1])*dy;
    V2z = (p3[2]-p1[2])*dz;

    Ax = V1y*V2z-V1z*V2y;
    Ay = V1z*V2x-V1x*V2z;
    Az = V1x*V2y-V1y*V2x;

    (*area) += 0.5*sqrt(Ax*Ax+Ay*Ay+Az*Az);

  }

  (*npolys) += num_o_polys;
}


/**************************** get_max_min ****************************/

void FTOC(get_max_min)(data,xdim,ydim,zdim,max,min)
register double *data;
int *xdim,*ydim,*zdim;
double *max,*min;
{

  double *enddata;

  enddata = data + ((*xdim) * (*ydim) * (*zdim));
  *max = *min = *(data++);
  for ( ;data<enddata;data++) {
    if (*data > *max)
       *max = *data;
    else
       if (*data < *min)
          *min = *data;
  }

}
