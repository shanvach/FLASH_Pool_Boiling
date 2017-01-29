Module VParticles_DPD_interface 
  implicit none
  
#include "Flash.h"
#include "constants.h" 
#include "Flash_mpi.h"

  interface 
     subroutine IO_writeDpd(particles,p_count,firstfileflag)
       !----------- Arguments list -------------------
       integer,INTENT(INOUT):: p_count
       real,dimension(NPART_PROPS,p_count),INTENT(INOUT)::particles
       logical,optional, INTENT(INOUT) :: firstfileflag
     end subroutine IO_writeDpd
    end interface
  interface 
    subroutine IO_writeDpd_parallel(particles,p_count,firstfileflag)	
!----------- Arguments list ------------------
  integer,INTENT(INOUT):: p_count
  real,dimension(NPART_PROPS,p_count),INTENT(INOUT)::particles
  logical,optional, INTENT(INOUT) :: firstfileflag
     end subroutine IO_writeDpd_parallel
  end interface     
    
  end Module VParticles_DPD_interface
