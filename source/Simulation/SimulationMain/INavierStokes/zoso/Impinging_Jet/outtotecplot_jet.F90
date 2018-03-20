! Subroutine outtotecplot
!
! Subroutine to write out to Tecplot data in binary form.
!
! ---------------------------------------------------------------------------

  subroutine outtotecplot_jet(mype,time,dt,istep,count, &
                          timer,blockList,blockCount,firstfileflag)


  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkPtr, &
    Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_getBlkBoundBox, &
    Grid_getBlkCenterCoords

  use ins_interface, only : ins_velgradtensor

  implicit none

#include "constants.h"
#include "Flash.h"
  include "Flash_mpi.h"


  integer, intent(in) :: mype,istep,count,&
                         firstfileflag
  integer, intent(in) :: blockCount
  integer, intent(in) :: blockList(MAXBLOCKS)
  real, intent(in)    :: time,dt,timer
  

  ! Local variables    
  integer :: numblocks,var,i,j,k,lb,nxc,nyc,nzc,ii,jj,kk
  character(27) :: filename
  character(6) :: index_lb,index_mype

  real xedge(NXB+2),xcell(NXB+2)
  real yedge(NYB+2),ycell(NYB+2)
  real zedge(NZB+2),zcell(NZB+2)
  real intsx(NXB+2), intsy(NYB+2), intsz(NZB+2)

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  real facevarxx(NXB+2*NGUARD+1,NYB+2*NGUARD,NZB+2*NGUARD),&
       facevaryy(NXB+2*NGUARD,NYB+2*NGUARD+1,NZB+2*NGUARD),&
       facevarzz(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD+1)

  real, dimension(NXB+1,NYB+1,NZB+1) :: tpu,tpv,tpw,tpp,&
            tpdudxcorn, tpdudycorn, tpdudzcorn,& 
            tpdvdxcorn, tpdvdycorn, tpdvdzcorn,&
            tpdwdxcorn, tpdwdycorn, tpdwdzcorn,&
            vortx,vorty,vortz,omg,             &
            Sxy,Syz,Sxz,Oxy,Oyz,Oxz,Qcr,divpp,TVtpp,&
            tempp,mdotp,tnlp,tnvp,nxp,nyp,nzp,tprds


  real*4 arraylb(NXB+2,NZB+2,1)

  real, dimension(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD) :: tpdudxc,&
        tpdudyc,tpdudzc,tpdvdxc,tpdvdyc,tpdvdzc,tpdwdxc,&
        tpdwdyc,tpdwdzc


  integer blockID

  real del(3),dx,dy,dz
  real, dimension(MDIM)  :: coord,bsize
  real ::  boundBox(2,MDIM)

  integer*4 TecIni,TecDat,TecZne,TecNod,TecFil,TecEnd
  integer*4 VIsdouble
  integer*4 Debug,ijk,Npts,NElm
  character*1 NULLCHR

!-----------------------------------------------------------------------
!                                                         TecPlot set-up
!-----------------------------------------------------------------------
  Debug     = 0
  VIsdouble = 0
  NULLCHR   = CHAR(0)

  ijk = (NXB+2)*(NZB+2)
!-----------------------------------------------------------------------

! -- filetime.XX --

  write(filename, '("IOData/data_time.", i4.4)') mype

  ! create/clear filetime.XX if time = 0
  if(firstfileflag .eq. 0) then
     open(unit=33, file=filename, status='replace')

#ifdef FLASH_GRID_UG
     write(33,*) 'NXB, NYB, NZB'
     write(33,'(3i4.1)') NXB, NYB, NZB    
#else
     write(33,*) 'NXB, NYB, NZB, interp. order (prolong, restrict)'
     write(33,'(3i4.1)') NXB, NYB, NZB
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

  write(filename,'("./IOData/temp.",i4.4,".",i4.4,".plt")') &
        count, mype

  i = TecIni('AMR3D'//NULLCHR,                            &
             'x y t'//NULLCHR,  &
             filename//NULLCHR,                           &
             './IOData/'//NULLCHR,                  &
             Debug,VIsdouble)

  intsx    = (/ (real(i), i=-1,NXB) /)
  intsy    = (/ (real(i), i=-1,NYB) /)
  intsz    = (/ (real(i), i=-1,NZB) /)

  call int2char(mype,index_mype)

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
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
     xcell = xedge + dx/2.0;

     yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
     ycell = yedge + dy/2.0;

     zedge = coord(KAXIS) - bsize(KAXIS)/2.0 + dz*intsz;
     zcell = zedge + dz/2.0;

     if(abs(ycell(2)-0.5*dy) .lt. 1e-15) then

     ! Write Block Results into data file:
     call int2char(lb,index_lb)

     i = TecZne(                                                       &
                'ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR,  &
                 NXB+2,NZB+2,1,                  &
                 'BLOCK'//NULLCHR,                                     &
                 CHAR(0))

     ! Write x:
     do j=1,NZB+2
        do i=1,NXB+2
           arraylb(i,j,1) = sngl(xcell(i))
        enddo
     enddo
     i = TecDat(ijk,arraylb,0)

     ! Write z:
     do j=1,NZB+2
        do i=1,NXB+2
           arraylb(i,j,1) = sngl(zcell(j))
        enddo
     enddo
     i = TecDat(ijk,arraylb,0)
     
     do j=1,NZB+2
        do i=1,NXB+2
        arraylb(i,j,1) = sngl(solnData(TEMP_VAR,NGUARD-1+i,1,NGUARD-1+j))
        end do
     end do
     i = TecDat(ijk,arraylb,0)

     end if

  enddo

  i = TecEnd()

  write(filename,'("./IOData/uvel.",i4.4,".",i4.4,".plt")') &
        count, mype

  i = TecIni('AMR3D'//NULLCHR,                            &
             'x y u'//NULLCHR,  &
             filename//NULLCHR,                           &
             './IOData/'//NULLCHR,                  &
             Debug,VIsdouble)

  intsx    = (/ (real(i), i=-1,NXB) /)
  intsy    = (/ (real(i), i=-1,NYB) /)
  intsz    = (/ (real(i), i=-1,NZB) /)

  call int2char(mype,index_mype)

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
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
     xcell = xedge + dx/2.0;

     yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
     ycell = yedge + dy/2.0;

     zedge = coord(KAXIS) - bsize(KAXIS)/2.0 + dz*intsz;
     zcell = zedge + dz/2.0;

     if(abs(ycell(2)-0.5*dy) .lt. 1e-15) then

     ! Write Block Results into data file:
     call int2char(lb,index_lb)

     i = TecZne(                                                       &
                'ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR,  &
                 NXB+2,NZB+2,1,                  &
                 'BLOCK'//NULLCHR,                                     &
                 CHAR(0))
    
     ! Write x:
     do j=1,NZB+2
        do i=1,NXB+2
           arraylb(i,j,1) = sngl(xedge(i))
        enddo
     enddo
     i = TecDat(ijk,arraylb,0)

     ! Write z:
     do j=1,NZB+2
        do i=1,NXB+2
           arraylb(i,j,1) = sngl(zcell(j))
        enddo
     enddo
     i = TecDat(ijk,arraylb,0)
     
     do j=1,NZB+2
        do i=1,NXB+2
        arraylb(i,j,1) = sngl(facexData(VELC_FACE_VAR,NGUARD-1+i,1,NGUARD-1+j))
        end do
     end do
     i = TecDat(ijk,arraylb,0)

     end if

  enddo

  i = TecEnd()

  write(filename,'("./IOData/vvel.",i4.4,".",i4.4,".plt")') &
        count, mype

  i = TecIni('AMR3D'//NULLCHR,                            &
             'x y v'//NULLCHR,  &
             filename//NULLCHR,                           &
             './IOData/'//NULLCHR,                  &
             Debug,VIsdouble)

  intsx    = (/ (real(i), i=-1,NXB) /)
  intsy    = (/ (real(i), i=-1,NYB) /)
  intsz    = (/ (real(i), i=-1,NZB) /)

  call int2char(mype,index_mype)

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
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
     xcell = xedge + dx/2.0;

     yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
     ycell = yedge + dy/2.0;

     zedge = coord(KAXIS) - bsize(KAXIS)/2.0 + dz*intsz;
     zcell = zedge + dz/2.0;

     if(abs(ycell(2)-0.5*dy) .lt. 1e-15) then

     ! Write Block Results into data file:
     call int2char(lb,index_lb)

     i = TecZne(                                                       &
                'ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR,  &
                 NXB+2,NZB+2,1,                  &
                 'BLOCK'//NULLCHR,                                     &
                 CHAR(0))

     ! Write x:
     do j=1,NZB+2
        do i=1,NXB+2
           arraylb(i,j,1) = sngl(xcell(i))
        enddo
     enddo
     i = TecDat(ijk,arraylb,0)

     ! Write z:
     do j=1,NZB+2
        do i=1,NXB+2
           arraylb(i,j,1) = sngl(zedge(j))
        enddo
     enddo
     i = TecDat(ijk,arraylb,0)
     
     do j=1,NZB+2
        do i=1,NXB+2
        arraylb(i,j,1) = sngl(facezData(VELC_FACE_VAR,NGUARD-1+i,1,NGUARD-1+j))
        end do
     end do
     i = TecDat(ijk,arraylb,0)

     end if

  enddo

  i = TecEnd()


  if (mype .eq. 0) then
  write(*,*) ''
  write(filename,'("./IOData/data.",i4.4,".**.plt")') &
        count
  write(*,*) '*** Wrote plotfile to ',filename,' ****'
  endif

  return

66    format(i4.4,g23.15,g23.15,i8.1,i5.1,g23.15)

  End subroutine outtotecplot_jet
