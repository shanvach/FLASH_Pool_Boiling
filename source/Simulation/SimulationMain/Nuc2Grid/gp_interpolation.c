#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "mangle_names.h"

#define SMALL 1.0E-10

void dposvx_(char *fact, char *uplo, int *n, int *nrhs, double *a, int *lda,
             double *af, int *ldaf, char *equed, double *s, double *b, int *ldb,
             double *x, int *ldx, double *rcond, double *ferr, double *berr,
             double *work, int *iwork, int *info);

static double *a, *x, *b, *af, *s, *work, *berr, *ferr, rcond;
static int *iwork, info;
static int Nparticles, Ndim, Ncells, Nproperties, Nrhs;
static double *Particle_Locations, *Particle_Properties, *Cell_Locations;
static double Sigma, *Mean, *Amplitude;
static double *Interpolant;

void allocate_arrays(void);
void fill_arrays(void);
void solve_linear_problems(void);
void calculate_gp_parameters(void);
void calculate_interpolants(void);
void free_arrays(void);

/***********************************************************************************

Furnish a Gaussian Process interpolation of properties recorded at given particle
location to given cell locations.

INPUT PARAMETERS:

*nparticles is the number of particles bearing data to be interpolated.
*ncells is the number of cells to which the interpolations are to be made.
*ndim is the number of physical dimensions.
*nproperties is the number of properties to be interpolated.

For 0 <= iparticles < *nparticles ; 0 <= icell < *ncells ;
    0 <= idim < *ndim  ;   0 <= iproperty < *nproperties : 

particle_locations[idim + iparticle * (*ndim)]
   is the idim-coordinate of the iparticle-th particle

particle_properties[iproperty + iparticle * (*nproperties)]
   is the value of the iproperty-th property of the iparticle-th particle

cell_locations[idim + icell * (*ndim)]
   is the idim-coordinate of the icell-th cell

sigma
  is an input parameter, containing the physical distance scale
  (sigma parameter) that is to be used in the squared-exponential
  covariance function.  Must be positive.  Probably should be set to
  some scale length suggested by the local hydrodynamic flow, like
  rho/|grad rho|.


OUTPUT PARAMETERS:

interpolant[ i + iproperty * 2 + icell * 2 * (*nproperties)]
   for i = 0: the interpolated value of the iproperty-th property to
              the icell-th cell                               
   for i = 1: the error estimate for this interpolant


mean[iproperty]
 contains the constant mean value of the Gaussian process model
 estimated from the particle data for the  iproperty-th property.
 

amplitude[iproperty]
  contains the amplitude parameter (\Sigma^2) used in the
  squared-exponential covariance function.  It is estimated from the
  particle data for the  iproperty-th property.

If called from a Fortran program,

    double precision particle_locations(ndim, nparticles)
    double precision particle_properties(nproperties, nparticles)
    double precision cell_locations(ndim, ncells)
    double precision sigma, mean(nproperties), amplitude(nproperties)
    double precision interpolant(2, nproperties, ncells)

    call gp_interpolation_(nparticles, ncells, ndim, nproperties,  &
                           particle_locations, particle_properties, &
                           cell_locations, sigma, mean, amplitude, interpolant)

Note the underscore after the subroutine name.  Note also the assumed dimensions of
the arrays, which must be allocated by the calling program to the exact dimensions.

The return value is the value of the parameter INFO returned by the LAPACK routine
dposvx.  It is 0 for success, > 0 for "something bad happened" (see the dposvx man
page).

***********************************************************************************/

int FTOC(gp_interpolation)(int *nparticles, int *ncells, int *ndim, int *nproperties,
                      double particle_locations[(*ndim) * (*nparticles)],
                      double particle_properties[(*nproperties) * (*nparticles)],
                      double cell_locations[(*ndim) * (*ncells)],
                      double *sigma, double *mean, double *amplitude, 
                      double interpolant[2 * (*nproperties) * (*ncells)])
{
    Nparticles = *nparticles;
    Ncells = *ncells;
    Ndim = *ndim;
    Nproperties = *nproperties;
    Nrhs = Ncells + Nproperties + 1;
    Particle_Locations = particle_locations;
    Particle_Properties = particle_properties;
    Cell_Locations = cell_locations;
    Sigma = *sigma;
    Mean = mean;
    Amplitude = amplitude;
    Interpolant = interpolant;

    allocate_arrays();
    fill_arrays();
    solve_linear_problems();
    assert( info >= 0 );
    if ( info > 0 ) return info;

    calculate_gp_parameters();
    calculate_interpolants();
    free_arrays();
    return 0;
}

/*---------------------------------------------------------------------------*/
/*                                                                           */
/* Allocate the arrays required for the solution of the SPD linear problem.  */
/*                                                                           */
/*---------------------------------------------------------------------------*/
void allocate_arrays(void)
{
    a = (double *) calloc(Nparticles*Nparticles, sizeof (double));
    af = (double *) calloc(Nparticles*Nparticles, sizeof (double));
    s = (double *) calloc(Nparticles, sizeof (double));
    x = (double *) calloc(Nparticles*Nrhs, sizeof (double));
    b = (double *) calloc(Nparticles*Nrhs, sizeof (double));
    ferr = (double *) calloc(Nrhs, sizeof (double));
    berr = (double *) calloc(Nrhs, sizeof (double));
    work = (double *) calloc(3*Nparticles, sizeof (double));
    iwork = (int *) calloc(Nparticles, sizeof(int));
}


/*---------------------------------------------------------------------------*/
/*                                                                           */
/* Free the temporary arrays.                                                */
/*                                                                           */
/*---------------------------------------------------------------------------*/
void free_arrays(void)
{
    free(a);
    free(af);
    free(s);
    free(x);
    free(b);
    free(ferr);
    free(berr);
    free(work);
    free(iwork);
}


/*---------------------------------------------------------------------------*/
/*                                                                           */
/* Fill in the required arrays.                                              */
/*                                                                           */
/*---------------------------------------------------------------------------*/
void fill_arrays(void)
{
    int i1, i2, idim, icell, iprop;
    double x1, x2, rsq, sigsq = Sigma * Sigma;

    for (i1 = 0 ; i1 < Nparticles ; i1++) {

        for (i2 = 0 ; i2 < Nparticles ; i2++) {
            /* Covariance matrix */
            for (idim = 0, rsq = 0.0 ; idim < Ndim ; idim++) {
                x1 = Particle_Locations[idim + i1 * Ndim];
                x2 = Particle_Locations[idim + i2 * Ndim];
                rsq += (x1 - x2) * (x1 - x2);
            }
            a[i2 + i1 * Nparticles] = exp(-rsq / sigsq);
        }
        a[i1 + i1 * Nparticles] += SMALL;

        /* Prediction coupling vectors */
        for ( icell = 0 ; icell < Ncells ; icell++) {

            for (idim = 0, rsq = 0.0 ; idim < Ndim ; idim++) {
                x1 = Particle_Locations[idim + i1 * Ndim];
                x2 = Cell_Locations[idim + icell * Ndim];
                rsq += (x1 - x2) * (x1 - x2);
            }
            b[i1 + Nparticles * icell] = exp(-rsq/sigsq);
        }

        /* Data */
        for ( iprop = 0 ; iprop < Nproperties ; iprop++ ) {

            b[i1 + Nparticles * (Ncells + iprop)] = 
                                Particle_Properties[iprop + i1 * Nproperties];
        }

        /*  "One" vector */
        b[i1 + Nparticles * (Ncells + Nproperties)] = 1.0;
    }
}

/*---------------------------------------------------------------------------*/
/*                                                                           */
/* Call Lapack to solve the required SPD linear problems.                    */
/*                                                                           */
/*---------------------------------------------------------------------------*/
void solve_linear_problems(void)
{
    char equed[] = {'N',0,0,0}, *fact = "N", *uplo = "U";

    FTOC(dposvx)(fact, uplo, &Nparticles, &Nrhs, a, &Nparticles, af,  &Nparticles,
            equed, s, b, &Nparticles, x, &Nparticles, &rcond, ferr, berr,
            work, iwork, &info);
}

/*---------------------------------------------------------------------------*/
/*                                                                           */
/* Calculate the mean and amplitude GP parameters from the data.             */
/*                                                                           */
/*---------------------------------------------------------------------------*/
void calculate_gp_parameters(void)
{
    int ipart, iprop;
    double u_qinv_u, u_qinv_f, f_qinv_f;

    for ( iprop = 0 ; iprop < Nproperties ; iprop++) {
        u_qinv_u = u_qinv_f = f_qinv_f = 0.0;
        for ( ipart = 0 ; ipart < Nparticles ; ipart++) {

           u_qinv_u += x[ipart + Nparticles * (Ncells + Nproperties)];
           u_qinv_f += x[ipart + Nparticles * (Ncells + iprop)];
           f_qinv_f += x[ipart + Nparticles * (Ncells + iprop)]
                        * b[ipart + Nparticles * (Ncells + iprop)];
        }
        Mean[iprop] = u_qinv_f / u_qinv_u;
        Amplitude[iprop] = (f_qinv_f - u_qinv_f * u_qinv_f / u_qinv_u)
                            / ((double) Nparticles); 
    }
}

/*---------------------------------------------------------------------------*/
/*                                                                           */
/* Calculate the values of the interpolants and their errors.                */
/*                                                                           */
/*---------------------------------------------------------------------------*/
void calculate_interpolants(void)
{
    int icell, iprop, ipart;
    double q_qinv_q, q_qinv_f;

    for ( icell = 0 ; icell < Ncells ; icell++ ) {

        for ( ipart = 0, q_qinv_q = 0.0 ; ipart < Nparticles ; ipart++)    
            q_qinv_q += x[ipart + Nparticles * icell] *
                        b[ipart + Nparticles * icell];
                                             
        for ( iprop = 0 ; iprop < Nproperties ; iprop++ ) {
        
            for ( ipart = 0, q_qinv_f = 0.0 ; ipart < Nparticles ; ipart++)
                q_qinv_f += x[ipart + Nparticles * icell] *
                            ( b[ipart + Nparticles * (Ncells + iprop)] -
                              Mean[iprop] );

            Interpolant[0 + iprop * 2 + icell * 2*Nproperties] =
                                                    Mean[iprop] + q_qinv_f;

            Interpolant[1 + iprop * 2 + icell * 2*Nproperties] = sqrt(
                                                    Amplitude[iprop] *
                                                    ( 1.0 - q_qinv_q ) );
        }
    }
}
