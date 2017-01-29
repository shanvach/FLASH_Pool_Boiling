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

int gp_interpolation(int *nparticles, int *ncells, int *ndim, int *nproperties,
                      double particle_locations[(*ndim) * (*nparticles)],
                      double particle_properties[(*nproperties) * (*nparticles)],
                      double cell_locations[(*ndim) * (*ncells)],
                      double *sigma, double *mean, double *amplitude, 
                      double interpolant[2 * (*nproperties) * (*ncells)]);
