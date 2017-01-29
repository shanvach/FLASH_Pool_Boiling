! Subroutine outtotecplot
!
! Subroutine to write out to Tecplot data in binary form.
!
! ---------------------------------------------------------------------------
#include "constants.h"
#include "Flash.h"

  subroutine outtotecplot_uv(mype,time,dt,istep,count, &
                          timer,blockList,blockCount,firstfileflag)


  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkPtr, &
                 Grid_releaseBlkPtr, Grid_getBlkIndexLimits, &
                 Grid_getBlkBoundBox,Grid_getBlkCenterCoords

  use ins_interface, only : ins_velgradtensor

#ifdef FLASH_GRID_PARAMESH
  use physicaldata, only : interp_mask_unk, interp_mask_unk_res
#endif
  implicit none
#include "Flash_mpi.h"
  integer, intent(in) :: mype,istep,count,&
                         firstfileflag
  integer, intent(in) :: blockCount
  integer, intent(in) :: blockList(MAXBLOCKS)
  real, intent(in)    :: time,dt,timer
  

  ! Local variables    
  integer :: numblocks,var,i,j,k,lb,nxc,nyc,nzc
  character(27) :: filename
  character(6) :: index_lb,index_mype

  real xedge(NXB+1),xcell(NXB+1)
  real yedge(NYB+1),ycell(NYB+1)
  real zedge(NZB+1),zcell(NZB+1)
  real intsx(NXB+1), intsy(NYB+1), intsz(NZB+1)

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  real, dimension(NXB,NYB,NZB) :: tpp
  real, dimension(NXB+1,NYB,NZB) :: tpu
  real, dimension(NXB,NYB+1,NZB) :: tpv
  real, dimension(NXB,NYB,NZB+1) :: tpw

  real*4 arraylbp(NXB,NYB,NZB)
  real*4 arraylbu(NXB+1,NYB,NZB)
  real*4 arraylbv(NXB,NYB+1,NZB)
  real*4 arraylbw(NXB,NYB,NZB+1)
 
  integer blockID

  real del(3),dx,dy,dz
  real, dimension(MDIM)  :: coord,bsize
  real ::  boundBox(2,MDIM)

  integer*4 TecIni,TecDat,TecZne,TecNod,TecFil,TecEnd
  integer*4 VIsdouble
  integer*4 Debug,ijk,Npts,NElm,ijku,ijkv,ijkw,ijkp
  character*1 NULLCHR

  logical :: file_exists

!-----------------------------------------------------------------------
!                                                         TecPlot set-up
!-----------------------------------------------------------------------
  Debug     = 0
  VIsdouble = 0
  NULLCHR   = CHAR(0)
  ijkp       = (NXB)*(NYB)*(NZB)
  ijku       = (NXB+1)*(NYB)*(NZB)
  ijkv       = (NXB)*(NYB+1)*(NZB)
  ijkw       = (NXB)*(NYB)*(NZB+1)
!-----------------------------------------------------------------------


! -- filetime.XX --

  write(filename, '("IOData/dauv_time.", i4.4)') mype
  INQUIRE(FILE=filename, EXIST=file_exists)

  ! create/clear filetime.XX if time = 0
  if( (firstfileflag .eq. 0) .or. (.not. file_exists) ) then
     open(unit=33, file=filename, status='unknown')

#ifdef FLASH_GRID_UG
     write(33,*) 'NXB, NYB, NZB'
     write(33,'(3i4.1)') NXB, NYB, NZB    
#else
     write(33,*) 'NXB, NYB, NZB, interp. order (prolong, restrict)'
     write(33,'(5i4.1)') NXB, NYB, NZB, interp_mask_unk(1), &
             interp_mask_unk_res(1)         
#endif

     write(33,'(a23,a43,a49,a12)') 'file number, time, dt, ',&
             'step number, ',&
             'total number of blocks, number of blocks output, ',&
             'elapsed time'
         
     close(33)
  endif

 ! write timestep data to filetime.XX on each processor

  open(unit=33, file=filename, status='old', position='append')
  write(33,66) count, time, dt, istep,blockcount,timer
  close(33)


  ! -- data.XXXX.XX --
  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1
  nzc = NZB + NGUARD + 1

! write UEVL to UEVL.XXXX.XX

  write(filename,'("./IOData/UVEL.",i4.4,".",i4.4,".plt")') &
        count, mype


  i = TecIni('AMR3D'//NULLCHR,                            &
             'x y z u'//NULLCHR,  &
             filename//NULLCHR,                           &
             './IOData/'//NULLCHR,                  &
             Debug,VIsdouble)

  intsx    = (/ (real(i), i=0,NXB) /)
  intsy    = (/ (real(i), i=0,NYB) /)
  intsz    = (/ (real(i), i=0,NZB) /)

  write(index_mype,"(I6.6)") mype

  do lb = 1,blockcount
  
     blockID =  blockList(lb)      

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)
     dx = del(IAXIS)
     dy = del(JAXIS)
     dz = del(KAXIS)

     ! Get Coord and Bsize for the block:
     ! Bounding box:
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)

     call Grid_getBlkCenterCoords(blockId,coord)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,facexData,FACEX)

     xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
     xcell = xedge(:) + dx/2.0;
    
     yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
     ycell = yedge(:) + dy/2.0;
    
     zedge = coord(KAXIS) - bsize(KAXIS)/2.0 + dz*intsz;
     zcell = zedge(:) + dz/2.0; 


     tpu = facexData(VELO_FACE_VAR, &
             NGUARD+1:NGUARD+NXB+1,NGUARD+1:NGUARD+NYB,NGUARD+1:NGUARD+NZB)
 
     ! Write Block Results into data file:
     write(index_lb,"(I6.6)") blockID
     i = TecZne(                                                       &
                'ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR,  &
                 NXB+1,NYB,NZB,                                    &
                 'BLOCK'//NULLCHR,                                     &
                 CHAR(0))

            
     ! Write x:
     do k=1,NZB
        do j=1,NYB
           do i=1,NXB+1
              arraylbu(i,j,k) = real(xedge(i),KIND=4) !sngl(xedge(i))
           enddo
        enddo
     enddo
     i = TecDat(ijku,arraylbu,0)


     ! Write y:
     do k=1,NZB
        do j=1,NYB
           do i=1,NXB+1
              arraylbu(i,j,k) = real(ycell(j),KIND=4) !sngl(yedge(j))
           enddo
        enddo
     enddo
     i = TecDat(ijku,arraylbu,0)


     ! Write z:
     do k=1,NZB
        do j=1,NYB
           do i=1,NXB+1
              arraylbu(i,j,k) = real(zcell(k),KIND=4) !sngl(zedge(k))
           enddo
        enddo
     enddo
     i = TecDat(ijku,arraylbu,0)


     ! Write u:
     arraylbu = real(tpu,KIND=4) !sngl(tpu)
     i = TecDat(ijku,arraylbu,0)

  enddo

  i = TecEnd()
  call Grid_releaseBlkPtr(blockID,facexData,FACEX)

  if (mype .eq. 0) then
  write(*,*) ''
  write(filename,'("./IOData/UVEL.",i4.4,".**.plt")') &
        count
  write(*,*) '*** Wrote plotfile to ',filename,' ****'
  endif

! write VEVL to VEVL.XXXX.XX

  write(filename,'("./IOData/VVEL.",i4.4,".",i4.4,".plt")') &
        count, mype


  i = TecIni('AMR3D'//NULLCHR,                            &
             'x y z v'//NULLCHR,  &
             filename//NULLCHR,                           &
             './IOData/'//NULLCHR,                  &
             Debug,VIsdouble)

  intsx    = (/ (real(i), i=0,NXB) /)
  intsy    = (/ (real(i), i=0,NYB) /)
  intsz    = (/ (real(i), i=0,NZB) /)

  write(index_mype,"(I6.6)") mype

  do lb = 1,blockcount
  
     blockID =  blockList(lb)      

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)
     dx = del(IAXIS)
     dy = del(JAXIS)
     dz = del(KAXIS)

     ! Get Coord and Bsize for the block:
     ! Bounding box:
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)

     call Grid_getBlkCenterCoords(blockId,coord)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

     xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
     xcell = xedge(:) + dx/2.0;
    
     yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
     ycell = yedge(:) + dy/2.0;
    
     zedge = coord(KAXIS) - bsize(KAXIS)/2.0 + dz*intsz;
     zcell = zedge(:) + dz/2.0; 


     tpv = faceyData(IBLK_FACE_VAR, &
             NGUARD+1:NGUARD+NXB,NGUARD+1:NGUARD+NYB+1,NGUARD+1:NGUARD+NZB)

 
     ! Write Block Results into data file:
     write(index_lb,"(I6.6)") blockID
     i = TecZne(                                                       &
                'ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR,  &
                 NXB,NYB+1,NZB,                                    &
                 'BLOCK'//NULLCHR,                                     &
                 CHAR(0))

            
     ! Write x:
     do k=1,NZB
        do j=1,NYB+1
           do i=1,NXB
              arraylbv(i,j,k) = real(xcell(i),KIND=4) !sngl(xedge(i))
           enddo
        enddo
     enddo
     i = TecDat(ijkv,arraylbv,0)


     ! Write y:
     do k=1,NZB
        do j=1,NYB+1
           do i=1,NXB
              arraylbv(i,j,k) = real(yedge(j),KIND=4) !sngl(yedge(j))
           enddo
        enddo
     enddo
     i = TecDat(ijkv,arraylbv,0)


     ! Write z:
     do k=1,NZB
        do j=1,NYB+1
           do i=1,NXB
              arraylbv(i,j,k) = real(zcell(k),KIND=4) !sngl(zedge(k))
           enddo
        enddo
     enddo
     i = TecDat(ijkv,arraylbv,0)


     ! Write v:
     arraylbv = real(tpv,KIND=4) !sngl(tpu)
     i = TecDat(ijkv,arraylbv,0)

  enddo

  i = TecEnd()
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

  if (mype .eq. 0) then
  write(*,*) ''
  write(filename,'("./IOData/VVEL.",i4.4,".**.plt")') &
        count
  write(*,*) '*** Wrote plotfile to ',filename,' ****'
  endif

! write WEVL to WEVL.XXXX.XX

  write(filename,'("./IOData/WVEL.",i4.4,".",i4.4,".plt")') &
        count, mype


  i = TecIni('AMR3D'//NULLCHR,                            &
             'x y z w'//NULLCHR,  &
             filename//NULLCHR,                           &
             './IOData/'//NULLCHR,                  &
             Debug,VIsdouble)

  intsx    = (/ (real(i), i=0,NXB) /)
  intsy    = (/ (real(i), i=0,NYB) /)
  intsz    = (/ (real(i), i=0,NZB) /)

  write(index_mype,"(I6.6)") mype

  do lb = 1,blockcount
  
     blockID =  blockList(lb)      

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)
     dx = del(IAXIS)
     dy = del(JAXIS)
     dz = del(KAXIS)

     ! Get Coord and Bsize for the block:
     ! Bounding box:
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)

     call Grid_getBlkCenterCoords(blockId,coord)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     tpw = 0.

     xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
     xcell = xedge(:) + dx/2.0;
    
     yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
     ycell = yedge(:) + dy/2.0;
    
     zedge = coord(KAXIS) - bsize(KAXIS)/2.0 + dz*intsz;
     zcell = zedge(:) + dz/2.0; 


     !tpw = facezData(VELO_FACE_VAR,:,:,:)
     tpw = facezData(VELC_FACE_VAR, &
            NGUARD+1:NGUARD+NXB, NGUARD+1:NGUARD+NYB, NGUARD+1:NGUARD+NZB+1)

 
     ! Write Block Results into data file:
     write(index_lb,"(I6.6)") blockID
     i = TecZne(                                                       &
                'ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR,  &
                 NXB,NYB,NZB+1,                                    &
                 'BLOCK'//NULLCHR,                                     &
                 CHAR(0))

            
     ! Write x:
     do k=1,NZB+1
        do j=1,NYB
           do i=1,NXB
              arraylbw(i,j,k) = real(xcell(i),KIND=4) !sngl(xedge(i))
           enddo
        enddo
     enddo
     i = TecDat(ijkw,arraylbw,0)


     ! Write y:
     do k=1,NZB+1
        do j=1,NYB
           do i=1,NXB
              arraylbw(i,j,k) = real(ycell(j),KIND=4) !sngl(yedge(j))
           enddo
        enddo
     enddo
     i = TecDat(ijkw,arraylbw,0)


     ! Write z:
     do k=1,NZB+1
        do j=1,NYB
           do i=1,NXB
              arraylbw(i,j,k) = real(zedge(k),KIND=4) !sngl(zedge(k))
           enddo
        enddo
     enddo
     i = TecDat(ijkw,arraylbw,0)


     ! Write u:
     arraylbw = real(tpw,KIND=4) !sngl(tpu)
     i = TecDat(ijkw,arraylbw,0)

  enddo

  i = TecEnd()
  call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  if (mype .eq. 0) then
  write(*,*) ''
  write(filename,'("./IOData/WVEL.",i4.4,".**.plt")') &
        count
  write(*,*) '*** Wrote plotfile to ',filename,' ****'
  endif

! write PRES to PRES.XXXX.XX

  write(filename,'("./IOData/PRES.",i4.4,".",i4.4,".plt")') &
        count, mype


  i = TecIni('AMR3D'//NULLCHR,                            &
             'x y z p'//NULLCHR,  &
             filename//NULLCHR,                           &
             './IOData/'//NULLCHR,                  &
             Debug,VIsdouble)

  intsx    = (/ (real(i), i=0,NXB) /)
  intsy    = (/ (real(i), i=0,NYB) /)
  intsz    = (/ (real(i), i=0,NZB) /)

  write(index_mype,"(I6.6)") mype

  do lb = 1,blockcount
  
     blockID =  blockList(lb)      

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)
     dx = del(IAXIS)
     dy = del(JAXIS)
     dz = del(KAXIS)

     ! Get Coord and Bsize for the block:
     ! Bounding box:
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)

     call Grid_getBlkCenterCoords(blockId,coord)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     tpp = 0.

     xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
     xcell = xedge(:) + dx/2.0;
    
     yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
     ycell = yedge(:) + dy/2.0;
    
     zedge = coord(KAXIS) - bsize(KAXIS)/2.0 + dz*intsz;
     zcell = zedge(:) + dz/2.0; 


     tpp = solnData(IBLK_VAR, &
            NGUARD+1:NGUARD+NXB, NGUARD+1:NGUARD+NYB, NGUARD+1:NGUARD+NZB)

 
     ! Write Block Results into data file:
     write(index_lb,"(I6.6)") blockID
     i = TecZne(                                                       &
                'ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR,  &
                 NXB,NYB,NZB,                                    &
                 'BLOCK'//NULLCHR,                                     &
                 CHAR(0))

            
     ! Write x:
     do k=1,NZB
        do j=1,NYB
           do i=1,NXB
              arraylbp(i,j,k) = real(xcell(i),KIND=4) !sngl(xedge(i))
           enddo
        enddo
     enddo
     i = TecDat(ijkp,arraylbp,0)


     ! Write y:
     do k=1,NZB
        do j=1,NYB
           do i=1,NXB
              arraylbp(i,j,k) = real(ycell(j),KIND=4) !sngl(yedge(j))
           enddo
        enddo
     enddo
     i = TecDat(ijkp,arraylbp,0)


     ! Write z:
     do k=1,NZB
        do j=1,NYB
           do i=1,NXB
              arraylbp(i,j,k) = real(zcell(k),KIND=4) !sngl(zedge(k))
           enddo
        enddo
     enddo
     i = TecDat(ijkp,arraylbp,0)


     ! Write u:
     arraylbp = real(tpp,KIND=4) !sngl(tpu)
     i = TecDat(ijkp,arraylbp,0)

  enddo

  i = TecEnd()
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  if (mype .eq. 0) then
  write(*,*) ''
  write(filename,'("./IOData/PRES.",i4.4,".**.plt")') &
        count
  write(*,*) '*** Wrote plotfile to ',filename,' ****'
  endif


66    format(i4.4,g23.15,g23.15,i8.1,i5.1,g23.15)

  return


  End subroutine outtotecplot_uv

