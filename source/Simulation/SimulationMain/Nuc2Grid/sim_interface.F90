Module sim_interface

  implicit none

  interface
     subroutine sim_readNucOutput(part_props,maxParticlesPerProc,numParticles,iPartFile)
       
       integer, intent(IN) :: part_props,maxParticlesPerProc,numParticles,iPartFile
       
     end subroutine sim_readNucOutput
  end interface


end Module sim_interface
