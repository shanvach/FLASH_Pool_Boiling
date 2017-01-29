!!****if* source/Grid/GridParticles/GridParticlesMove/UG/Directional/Grid_moveParticles
!!
!! NAME
!!  Grid_moveParticles
!!
!! SYNOPSIS
!!
!!  Grid_moveParticles(real(INOUT)    :: dataBuf(:,:),
!!                     integer(IN)    :: propCount 
!!                     integer(INOUT) :: localCount,
!!                     integer(IN)    :: maxCount,
!!                     integer(IN)    :: index_list(indexCount),
!!                     integer(IN)    :: indexCount
!!                     logical(IN)    :: coords_in_blk)
!!
!!  
!! DESCRIPTION 
!!  
!!  If some of the particles (active or tracer) need to move to the part of the physical
!!  domain that is mapped to a different processor, this function moves their 
!!  data structure to the appropriate processors
!!
!!  This implementations gets particles to the correct processor by first sweeping in the
!!  x direction and moving particles to the x direction neighbors.  Then another sweep 
!!  in the y direction moves particles to the y neighbors, followed by z.  A particle 
!!  then could physically be moved 3 times.  We proceeded in this way because this method
!!  reduces the number of exchanges which are expensive because of latency.
!!  
!! ARGUMENTS 
!!
!!  dataBuf : the data structure containing the data of interest
!!              It is two dimensional real array, the first dimension
!!              represents properties associated with the data 
!!              structure, and 
!!              second dimension is index to individual elements in 
!!              the datastructure.
!!
!! propCount : number of properties for this datastructure 
!!
!!  localCount : While coming in it contains the current 
!!               number of elements in the data structure mapped to
!!               this processor. After all the data structure 
!!               movement, the number might change, 
!!               and the new value is put back into it
!!  maxCount : This is parameter determined at runtime, 
!!             and is the maximum count of elements 
!!             that a simulation expects to have. 
!!             All the arrays  are allocated based on this number
!!
!!  index_list : The list of fields in the incoming dataBuf,
!!               the list are used to make the indices of tag, block
!!               processor and physical location of each individual
!!               data item to the GridParticles subunit so it can
!!               move the concerned data item appropriately
!!  indexCount : The count of fields included in index_list
!!
!!  coords_in_blk   : if true then re-gridding just happened and blocks have changed
!!
!!***


#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif


subroutine Grid_moveParticles(dataBuf, propCount, maxCount, localCount, &
     index_list, indexCount, coords_in_blk)
 
  use Grid_data, ONLY : gr_axisMe,gr_axisNumProcs,gr_axisComm,gr_imin,gr_jmin,gr_kmin,gr_domainBC, &
       gr_meshAcrossMe,gr_meshAcrossNumProcs,gr_useParticles, gr_meshNumProcs,gr_meshMe

  use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_getDeltas
  use gr_ptInterface, ONLY : gr_ptMoveOffProc, gr_ptSetIndices, gr_ptResetIndices, gr_ptMoveSieve
 
  use gr_ptData, ONLY : gr_ptPosx, gr_ptPosy, gr_ptPosz, gr_ptProc,gr_ptblk,gr_ptSourceBuf,gr_ptDestBuf
  
  use Particles_data , ONLY: pt_posInitialized
  use gr_ptVPData, only : gr_ptVPBufferFactor
  implicit none

#include "constants.h"
#include "Flash.h"


  integer,intent(IN) :: maxCount, propCount, indexCount
  integer,dimension(indexCount), intent(IN) :: index_list
  integer,intent(INOUT) :: localCount

  real, dimension(propCount, maxCount),intent(INOUT) :: dataBuf
  logical, intent(IN) :: coords_in_blk



  integer i, j
  integer,dimension(MDIM) :: pos 
  real, dimension(MDIM)   :: deltas
  integer :: lnegh,rnegh
  integer :: blockID=1
  real,dimension(LOW:HIGH,MDIM) :: bndBox
  logical ::lbdry,rbdry

  integer :: MyPE,MoveOffBlk_index
  integer :: numDest,m
  real :: rc 
 
  
  call gr_ptSetIndices(index_list,indexCount)
  call Grid_getBlkBoundBox(blockID,bndBox)
  call Grid_getDeltas(blockID,deltas)
  
  if (.not.pt_posInitialized) then
     
     dataBuf(gr_ptBlk,localCount+1:maxCount)=NONEXISTENT
     
     if (coords_in_blk) then
        numDest=0
        ! Initialize the blk property so every particle is checked
        
        dataBuf(gr_ptBlk,1:localCount)=CONSTANT_ONE
        gr_ptSourceBuf(:,1:localCount)=dataBuf(:,1:localCount)
        ! HE: Zero the local count, because we don't now the number of particles 
        ! that actually belong to the current processor
        m=0
        call gr_ptLocalMatch(dataBuf,m,propCount,maxCount,gr_ptSourceBuf,&
             localCount, gr_ptDestBuf,numDest)
        localCount=m
       
        dataBuf(gr_ptBlk,localCount+1:maxCount)=NONEXISTENT
        
        call gr_ptMoveSieve(dataBuf,localCount,propCount,&
             maxCount, numDest, gr_ptProc)
        write(*,*) 'Partilces placed in their destinations using gr_ptMoveSieve'
	pt_posInitialized=.true.
        return
     end if
     
  end if
    
  pos(IAXIS)=gr_ptPosx
  pos(JAXIS)=gr_ptPosy
  pos(KAXIS)=gr_ptPosz
  
  
#ifdef DEBUG_GRID
  print *,'In Grid_moveParticles, basic geometry is found'
#endif

  do i = 1,NDIM

     lnegh=gr_axisMe(i)-1             !! determine the left negh along i
     lbdry = .false.    
     if(lnegh<0) then                 !! find out if the block is on 
        lnegh=gr_axisNumProcs(i)-1    !! left boundary
        lbdry = .true.
     end if
     
     rnegh=gr_axisMe(i)+1    !! do the same for right neighbor 
     rbdry = .false.         !! boundary
     if(rnegh==gr_axisNumProcs(i)) then
        rnegh=0
        rbdry = .true.
     end if
     rc=deltas(i)*gr_ptVPBufferFactor;
        
     call gr_ptMoveOffProc(LOW,i,pos(i), propCount, maxCount,lbdry,&
                           lnegh,rnegh,bndBox(:,i),localCount,dataBuf,rc)
        
#ifdef DEBUG_GRID
     print *,'In Grid_moveParticles, LOW end moved'
#endif     
     call gr_ptMoveOffProc(HIGH,i,pos(i), propCount, maxCount,rbdry,&
                           rnegh,lnegh,bndBox(:,i),localCount,dataBuf,rc)
     
  end do
  
#ifdef DEBUG_GRID
  print *,'In Grid_moveParticles, all done'
#endif
  call gr_ptResetIndices(index_list,indexCount)

end subroutine Grid_moveParticles
