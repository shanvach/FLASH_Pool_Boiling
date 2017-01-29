!!****if* source/Simulation/SimulationMain/VParticles_DPD/pt_initPositions
!!
!! NAME
!!    pt_initPositions
!!
!! SYNOPSIS
!!
!!    call pt_initPositions(integer(in)  :: blockID,
!!                          logical(out) :: success)
!!
!! DESCRIPTION
!!
!!    Initializes particle locations for one block in the grid.
!!
!! ARGUMENTS
!!
!!  blockID:        local block ID containing particles to create
!!
!!  success:        returns .TRUE. if positions for all particles
!!                  that should be assigned to this block have been
!!                  successfully initialized.
!!
!!***

#define NPART_PROPS_IO  7
subroutine pt_initPositions (blockID,success)
  use IO_data, ONLY : io_outputSplitNum, &
       io_comm, io_splitNumBlks, io_splitParts, io_globalMe,&
       io_meshNumProcs, io_meshMe
  use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_moveParticles
  use Grid_data, ONLY : gr_meshComm
  use Simulation_data,ONLY : sim_initPos,pt_NumPart,pt_NumBodies, &
       pt_NumBTypes,Aij,Connect,MAXBTYPES,pt_BodyTypes,MAXCONN
  use Particles_data, ONLY:  pt_numLocal, particles,  & 
       pt_maxPerProc, pt_xmin, pt_ymin, pt_zmin,pt_xmax, &
       pt_ymax, pt_zmax,pt_posAttrib,pt_velNumAttrib,    &
       pt_velAttrib,pt_typeInfo, pt_meshMe,pt_meshNumProcs, pt_meshComm, &
       pt_indexList, pt_indexCount
  use Driver_interface,ONLY: Driver_abortFlash

#ifdef USE_IO_C_INTERFACE
  use iso_c_binding, ONLY : c_loc
  use io_c_interface, ONLY : io_xfer_cont_slab
#else
#define c_loc(x) x
#endif

  !----------------
  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "Particles.h"
#include "Flash_mpi.h"
#include "io_flash.h"  
 
!--- Argument List
  integer, INTENT(in) :: blockID
  logical,intent(OUT) :: success
  
! ------------- local variables ------------
  real    :: xpos, ypos, zpos, rpos, bxl, byl, bzl, bxu, byu, bzu
  real, dimension(CONSTANT_TWO,MDIM):: boundBox
  real,dimension(pt_NumPart) :: ref_xp,ref_yp,ref_zp
  real,dimension(pt_NumPart) :: BodyType,parent,beadtype,internal_index
  real :: beadparam(MAXBTYPES)
  integer :: bufsize
  
  character(len=*), parameter :: dsetName = "tracer particles", &
       localnp_str = "localnp"

  integer :: parts_given,p,int1,int2,int3,int4
  integer :: part_props=NPART_PROPS,mapType,nlks,nnodes
  integer :: k,i,j,jj
  integer :: total_pt_local,ierr
  integer,allocatable :: Links(:,:) ,conn(:)
  integer,dimension(6):: SysData
  integer :: istat,reLocalNumParticles
  integer :: memOffset, fileOffset
  integer, parameter :: libType = IO_FILE_HDF5
  integer, parameter :: xferType = IO_READ_XFER
  integer, parameter :: memType = IO_FLASH_DOUBLE
  integer, parameter :: partArrayDims = 2
  integer, parameter :: dsetLen = len_trim(dsetName)
  integer :: typeMatchedXfer, err, particlesPerBlkSize
  integer :: fileID
  integer :: localOffset, splitOffset, partSplitOffset, localPartOffset
  integer :: particleOffset

  logical :: IsInBlock, coords_in_blk
  real, allocatable :: particles2(:,:)
  
  !-----------------------------------------
  !write(*,*)'NPART_PROP_IO',NPART_PROPS_IO
  p = pt_numLocal;
  
  !--- The particles file are always double precision  
  typeMatchedXfer = 0
  
  !write(*,*)'Initial number of particles on this proc= ',p;
  
  ! Read the simulation and system data from MEMBRANE.dat file
  
  if (pt_meshMe==MASTER_PE) then 
     open (unit=9, file = 'MEMBRANE.dat', status = 'old', action = 'read')
     read (9, *)int1,int2,int3,int4
     if (int1/=pt_NumPart)   call Driver_abortFlash(" Check the dat file and make sure it has the correct number of particles defined in flash.par")
     if (int2/=pt_NumBTypes) call Driver_abortFlash(" Check the dat file and make sure it has the correct number of particles types defined in flash.par")
     if (int3/=pt_BodyTypes) call Driver_abortFlash(" Check the dat file and make sure it has the correct number of bodies defined in flash.par")  
     if (int4/=pt_NumBodies) call Driver_abortFlash(" Check the dat file and make sure it has the correct number of bodies defined in flash.par") 
     
     write(*,*)'Total number of partilces in data file =',pt_NumPart,  ' containing',pt_NumBTypes,'bead type'
     write(*,*)'Total number of bodies in the data file=',pt_NumBodies,' containing',pt_BodyTypes,'body types'

     write(*,*)pt_numbtypes
     do i=1,pt_numbtypes
        read (9, *)beadparam(1:pt_numbtypes)
        Aij(i,1:pt_numbtypes)=beadparam(1:pt_numbtypes)
     end do
     
     ! Reading connectivity topology
     do i=1,pt_BodyTypes
        read (9, *) nlks
        write(*,*)  'Number of links in bodytype',i,' is ',nlks
        allocate(Links(nlks,2))
        
        
        do j=1,nlks
           read(9,*) Links(j,1:2)
        end do
        
        nnodes=maxval(Links)
        write(*,*)'No. of beads in body:',i,'is: ',nnodes
        
     end do
     
     close(9)
     SysData=(/pt_NumPart,pt_numbtypes,pt_numBodies,pt_BodyTypes,nlks,nnodes/);
  end if
  
  
  
  ! Broadcast the Aij, nlks, numbtypes, Connect, BodyTypes
  call MPI_Bcast (SysData, 6, FLASH_INTEGER, MASTER_PE, &
       &                  pt_MeshComm, istat)
  
  if (pt_meshMe/=MASTER_PE) then
     pt_NumPart=SysData(1);
     pt_numbtypes=SysData(2);
     pt_numBodies=SysData(3);
     pt_BodyTypes=SysData(4);
     nlks=SysData(5);
     nnodes=SysData(6);
  end if
  
  if (.NOT.allocated(Connect)) allocate(Connect(pt_BodyTypes))
  if (.NOT.allocated(Links))   allocate(Links(nlks,2))
  
  bufsize=nlks*2;
  
  call MPI_Bcast (Links(:,:) , bufsize , FLASH_INTEGER, MASTER_PE, &
       &                  pt_MeshComm, istat)
!!$  call MPI_Bcast (Links(CONSTANT_ONE:nlks,CONSTANT_TWO), nlks, FLASH_INTEGER, MASTER_PE, &
!!$       &                  pt_MeshComm, istat)
  
  ! Fill the connectivity arrays.
  do i = 1, pt_BodyTypes
     Connect(i)%numLinks = nlks   
     
     if (.NOT.allocated(Connect(i)%links))allocate(Connect(i)%links(nnodes,MAXCONN))
     if (.NOT.allocated(conn))allocate(conn(nnodes))
     
     Connect(i)%links(:,:)=0;
     conn(:)=1;
     do jj=1,nlks
        Connect(i)%links(Links(jj,1),conn(Links(jj,1)))=Links(jj,2)
        conn(Links(jj,1))=conn(Links(jj,1))+1
     end do
     
  end do

  
  ! Call hdf5 reader to read the particles from the hdf5 files
  call io_h5open_file_for_read(FileID,'MEMBRANE.h5 ', io_comm, io_outputSplitNum)
  
  reLocalNumParticles=pt_NumPart/pt_meshNumProcs
  write(*,*)'reLocalNumParticles',reLocalNumParticles
  if(pt_meshMe==(pt_meshNumProcs-1))&
       reLocalNumParticles=pt_NumPart-reLocalNumParticles*pt_meshMe
  
  if (reLocalNumParticles > pt_maxPerProc) then
     call Driver_abortFlash &
          ('[io_ptReadParticleData] ERROR: too many particles on this proc; increase pt_maxPerProc')
  end if
  
  !allocate particles data structure
  if (.NOT.allocated(particles2)) &
       allocate (particles2(NPART_PROPS_IO,pt_maxPerProc), stat=ierr)
  if (ierr /= 0) &
       call Driver_abortFlash("io_ptReadParticleData:  could not allocate particle array")
  
  !particles must be initialized or the entire particles algorithm will fail
  particles = NONEXISTENT
  particles2 = NONEXISTENT
  !write(*,*) 'No of processors:',pt_meshNumProcs
  if (pt_meshNumProcs>1) then 
     
     !now get the particle offset
     call io_getParticleOffset( reLocalNumParticles, io_splitParts, particleOffset)
     
     call MPI_ALLREDUCE(particleOffset, partSplitOffset, 1, FLASH_INTEGER, &
          MPI_MIN, io_comm, ierr)
     localPartOffset = particleOffset - partSplitOffset
  else
     localPartOffset=0
  end if

  write(*,*) 'The localpartOffset=',localPartOffset 
  call io_xfer_cont_slab(io_globalMe, &
       FileID, &
       libType, &
       xferType, &
       typeMatchedXfer, &
       dsetName, &
       dsetLen, &
       memType, &
       (/pt_maxPerProc,NPART_PROPS_IO/), &
       (/0,0/), &
       (/reLocalNumParticles,NPART_PROPS_IO/), &
       (/localPartOffset,0/), &
       (/reLocalNumParticles,NPART_PROPS_IO/), &
       partArrayDims, &
       c_loc(particles2(1,1)), err)
  if (err /= 0)call Driver_abortFlash("Error reading particles")
  
  !CD: Now we need copy across the initialized locations into the
  !real particles data structure.
  particles(POSX_PART_PROP,1:reLocalNumParticles) = particles2(1,1:reLocalNumParticles)
  particles(POSY_PART_PROP,1:reLocalNumParticles) = particles2(2,1:reLocalNumParticles)
  particles(POSZ_PART_PROP,1:reLocalNumParticles) = particles2(3,1:reLocalNumParticles)
  !HE: other particles properties loaded to the particles data structure
  particles(PRNT_PART_PROP,1:reLocalNumParticles) = particles2(4,1:relocalNumParticles)
  particles( BDT_PART_PROP,1:reLocalNumParticles) = particles2(5,1:relocalNumParticles)
  particles(INTR_PART_PROP,1:reLocalNumParticles) = particles2(6,1:relocalNumParticles)
  particles(BDYT_PART_PROP,1:reLocalNumParticles) = particles2(7,1:relocalNumParticles)  
  particles(FNPX_PART_PROP:FNPZ_PART_PROP,:) = 0;
  particles(FNX_PART_PROP:FNZ_PART_PROP,:)   = 0;
  particles(VELX_PART_PROP:VELZ_PART_PROP,:) = 0;
  particles(VNPX_PART_PROP:VNPZ_PART_PROP,:) = 0;
  particles(VIX_PART_PROP:VIZ_PART_PROP,:)   = 0;
  particles(BLK_PART_PROP,:) = NONEXISTENT

 
!!$  do i=1,reLocalNumParticles
!!$     write(*,*) particles(POSX_PART_PROP:POSZ_PART_PROP,i)
!!$  end do
 
  !CD: Set pt_numLocal which is a particles unit variable keeping track
  !of the number of particles currently resident on this processor.
  pt_numLocal = reLocalNumParticles 

  
  ! HE: Just for reference (properties order in the input file MEMBRANE.h5) 
  ! read (9, *)ref_xp(i),ref_yp(i),ref_zp(i),parent(i),beadtype(i),internal_index(i),BodyType(i) 
  
  bufsize=MAXBTYPES* MAXBTYPES; 
  
  call MPI_Bcast (Aij(:,:), bufsize , FLASH_REAL, MASTER_PE, &
       &                  pt_meshComm, istat)
  !if (istat== MPI_SUCCESS) write(*,*)'Successful broadcasting'

  mapType=pt_typeInfo(PART_MAPMETHOD,1)
  
  if (allocated(conn))       deallocate(conn)
  if (allocated(particles2)) deallocate(particles2)
 
  !-----------
    
  coords_in_blk=.true.
  call Grid_moveParticles(particles,NPART_PROPS,pt_maxPerProc,&
       pt_numLocal,pt_indexList, pt_indexCount, coords_in_blk) 
  write(*,*) '  '
  write(*,*) '  '
  write(*,*) '  '
  write(*,*) '  '
  write(*,*) '  '
  write(*,*) '  '
  write(*,*) '  '
  write(*,*) '  '

!!$  do i=1,pt_numlocal
!!$     write(*,*) i,particles(BLK_PART_PROP,i)
!!$  end do
  print *, "Local num. particles after moving particles to the correct destination =", &
       pt_numLocal
  
  success = .TRUE.
  
  
  !----------------------------------------------------------------------
  
end subroutine pt_initPositions


