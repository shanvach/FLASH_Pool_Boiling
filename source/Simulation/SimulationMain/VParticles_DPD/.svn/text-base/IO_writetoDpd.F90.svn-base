subroutine IO_writeDpd(particles,p_count,firstfileflag)
  use HDF5
  use Driver_data, ONLY: dr_globalMe, dr_nbegin, &
       dr_nend, dr_dt, dr_wallClockTimeLimit, &
       dr_tmax, dr_simTime, dr_redshift, &
       dr_nstep, dr_dtOld, dr_dtNew, dr_nbegin, dr_restart
  use Simulation_data, ONLY : domainsize, pt_NumPart
  use Particles_data, ONLY: pt_maxPerProc
  use Grid_data,ONLY : gr_globalDomain, gr_meshNumProcs 
  use Grid_interface, ONLY:Grid_sortParticles

  implicit none

#include "Flash.h"
#include "constants.h" 


  !----------- Arguments list -------------------
  integer,INTENT(INOUT):: p_count
  real,dimension(NPART_PROPS,p_count),INTENT(INOUT)::particles
  logical,optional, INTENT(INOUT) :: firstfileflag
  !----------------------------------------------
  
  ! Local variables
  real,dimension(NDIM,p_count) ::pos,v,f 
  real :: tages(1:p_count)
  integer :: count,i
  character(len=6) :: index_count,index_proc
  character(len=4) :: name
  integer,dimension(MDIM):: opind,npind,ovind,intvind,nvind,ofind,nfind
  integer, dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk
  logical :: writexyz=.true.
  !----------------------------------------------
  

  if(.not.present(firstfileflag))firstfileflag=.false.

  ! Create the file name      
  call int2char(dr_nstep,index_count)
  call int2char(dr_globalMe,index_proc)
  open(unit=113,file='./IOData/part_'//index_count//'.'//index_proc//'.txt',form='formatted')
  
  call pt_dpdSetIndices(opind,npind,ovind,nvind,intvind,ofind,nfind)
  ! Read the data from the particles array to update the positions
  tages(:)  = particles(TAG_PART_PROP,1:p_count);
  do i=1,NDIM
     !f(i,:)      = particles(nfind(i),1:p_count);
     v(i,:)      = particles(nvind(i),1:p_count);
     pos(i,:)    = particles(opind(i),1:p_count);
     
  end do
  write(113,'(3(I7,1X),2F16.10)')pt_NumPart,p_count,int(dr_nstep),dr_dt,dr_Simtime
  do i=1,p_count 
#if NDIM==2

     write(113,'(I5,4F16.10)')int(tages(i)),pos(1:NDIM,i),v(1:NDIM,i)
#elif NDIM==3

     write(113,'(I7,2X,6(F16.10,1X))')int(tages(i)),pos(1:NDIM,i),v(1:NDIM,i)
#endif
 
  end do
  close(113);

  open(unit=115,file='./IOData/part'//index_count//'.'//index_proc//'.plt',form='formatted')
  
  write(115,'(A)') 'VARIABLES = "XPOS" , "YPOS", "ZPOS"'
  write(115,'(A,I8,A)')'ZONE I=',p_count,',DATAPACKING = POINT' 
  
  do i = 1,p_count
#if NDIM==3
     !write(6,'(3(2X,F16.12))') pos(1:3,i) 
     write(115,'(3(2X,F16.12))') pos(1:3,i)  !,pos(2,i),pos(3,i)
#endif
  enddo
  close(115)
  !stop
  !write(*,*) '2 File flag in IO_writeDpd', firstfileflag,'proc#',dr_globalMe

  !----------- Write the xyz file ---------------------
  if (writexyz) then
     if (gr_meshNumProcs==1) then
        if(firstfileflag) then 
           open(unit=117,file='./IOData/parts.xyz',form='formatted',STATUS='REPLACE')
           
        else
           open(unit=117,file='./IOData/parts.xyz',form='formatted',POSITION='APPEND')
        end if
        write(117,*)p_count 
        write(117,*)'RBC'     
        
        do i = 1,p_count
           if     (particles(BDT_PART_PROP,i)==1.) then
              name = 'N'
           elseif (particles(BDT_PART_PROP,i)==2.) then
              name = 'P'
           elseif (particles(BDT_PART_PROP,i)==3.) then
              name = 'C'
           elseif (particles(BDT_PART_PROP,i)==4.) then
              name = 'H'
           end if
           
#if NDIM==2
           if (particles(TAG_PART_PROP,i)> 0 ) write(117,'(A,2(2X,F16.12))') name,pos(1:NDIM,i)  
#elif NDIM==3
           if (particles(TAG_PART_PROP,i)> 0 ) write(117,'(A,3(2X,F16.12))') name,pos(1:3,i)  
#endif
        end do
        
     close(117)
  else
      open(unit=117,file='./IOData/parts.'//index_count//'.'//index_proc//'.xyz',form='formatted')
      write(117,*)p_count 
        write(117,*)'RBC'     
        
        do i = 1,p_count
           if     (particles(BDT_PART_PROP,i)==1.) then
              name = 'N'
           elseif (particles(BDT_PART_PROP,i)==2.) then
              name = 'P'
           elseif (particles(BDT_PART_PROP,i)==3.) then
              name = 'C'
           elseif (particles(BDT_PART_PROP,i)==4.) then
              name = 'H'
           end if
           
#if NDIM==2
           if (particles(TAG_PART_PROP,i)> 0 ) write(117,'(A,2(2X,F16.12))') name,pos(1:NDIM,i)  
#elif NDIM==3
           if (particles(TAG_PART_PROP,i)> 0 ) write(117,'(A,3(2X,F16.12))') name,pos(1:3,i)  
#endif
        end do
        
     close(117)
  end if
end if
!write(*,*) '3 File flag in IO_writeDpd', firstfileflag,'proc#',dr_globalMe
end subroutine IO_writeDpd




#ifdef DE
        if ((pos(IAXIS,i)>gr_globalDomain(HIGH,IAXIS)).or.(pos(IAXIS,i)<gr_globalDomain(LOW,IAXIS))) then
           write(*,*)'This particle is out of the x boundary'
           write(*,*)'pos=',particles(POSX_PART_PROP,i),particles(TAG_PART_PROP,i)
           !stop
        elseif ((pos(JAXIS,i)>gr_globalDomain(HIGH,JAXIS)).or.(pos(JAXIS,i)<gr_globalDomain(LOW,JAXIS))) then
           write(*,*)'This particle is out of the y boundary'
           write(*,*)'pos=',particles(POSY_PART_PROP,i),particles(TAG_PART_PROP,i)
           !stop
        elseif ((pos(KAXIS,i)>gr_globalDomain(HIGH,KAXIS)).or.(pos(KAXIS,i)<gr_globalDomain(LOW,KAXIS))) then
           write(*,*)'This particle is out of the z boundary'
           write(*,*)'pos=',particles(POSZ_PART_PROP,i),particles(TAG_PART_PROP,i)
           !stop
        end if
#endif
