#include "Particles.h"
#include "Flash.h"

module Particles_sinkData
  public ::  MAX_MSGS, maxsinks, n_empty
  logical, save :: RunningParticles = .false.
  integer, save :: n_empty
  integer, save :: MaxParticlesPerProc
  logical, save :: UseSinkParticles

  integer, parameter :: maxsinks = 2048
  integer, save, pointer, dimension(:)  :: NumParticlesPerBlock
  integer, parameter :: MAX_MSGS = 12
  integer, parameter :: nrep_pbc = 2

  ! It's probably not the best way to do it, but for MPI purposes
  ! let's do a separate send_buff and recv_buff for each particle
  ! property
  !real, dimension(MAX_MSGS), pointer :: send_buff, recv_buff
  real,  dimension(MAX_MSGS) :: send_buff, recv_buff

  !logical, save, allocatable, dimension(:)         :: is_empty
  integer, save :: ipx, ipy, ipz, ipvx, ipvy, ipvz, ipm, iptag
  integer, save :: ipblk, iplx, iply, iplz, ipmdot, ipt, ipsc
  integer, save :: ipdtold, ipcpu, iold_pmass, ipraccr, ipmgas
  integer, save :: ipbflx, ipbfly, ipbflz

  integer, save :: pt_maxSinksPerProc
  integer, parameter :: pt_sinkParticleProps = NPART_PROPS

  ! particles_local and particles_global refer to 
  ! sink particles - the local list and the global list
  real, save, allocatable, dimension(:,:) :: particles_local
  real, save, allocatable, dimension(:,:) :: particles_global
  
  integer, save :: local_tag_number
  integer, save :: localnp, localnpf
end module
