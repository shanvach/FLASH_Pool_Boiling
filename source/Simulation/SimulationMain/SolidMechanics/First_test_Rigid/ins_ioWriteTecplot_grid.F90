! Subroutine ins_ioWriteTecplot_grid
!
!
! ---------------------------------------------------------------------------
#include "constants.h"
#include "Flash.h"


  subroutine ins_ioWriteTecplot_grid(mype,time,dt,istep,count, &
                     timer,blockList,blockCount,firstfileflag)


  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkPtr,       &
  Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_getBlkBoundBox, &
  Grid_getBlkCenterCoords

#ifdef FLASH_GRID_PARAMESH
  use physicaldata, only : interp_mask_unk, interp_mask_unk_res
#endif

  implicit none
  include "Flash_mpi.h"
  integer, intent(in) :: mype,istep,count,&
                         firstfileflag
  integer, intent(in) :: blockCount
  integer, intent(in) :: blockList(MAXBLOCKS)
  real, intent(in)    :: time,dt,timer
  

  ! Local variables    
  integer :: numblocks,var,i,j,k,lb,nxc,nyc,nzc
  character(29) :: filename
  character(6) :: index_lb,index_mype

  real xedge(NXB+1),xcell(NXB+1)
  real yedge(NYB+1),ycell(NYB+1)
  real zedge(NZB+1),zcell(NZB+1)
  real intsx(NXB+1), intsy(NYB+1), intsz(NZB+1)

  real*4 arraylb(NXB+1,NYB+1,NZB+1)
 
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
  ijk       = (NXB+1)*(NYB+1)*(NZB+1)
!-----------------------------------------------------------------------


! -- filetime.XXXX --

  write(filename, '("IOData/Grid_time.", i4.4)') mype

  ! create/clear filetime.XX if time = 0
  if(firstfileflag .eq. 0) then
     open(unit=33, file=filename, status='replace')

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

! write solution data to Grid.XXXXXX.XXXX

  write(filename,'("./IOData/Grid.",i6.6,".",i4.4,".plt")') &
        count, mype


  i = TecIni('AMR3D'//NULLCHR,                            &
             'x y z'//NULLCHR,  &
             filename//NULLCHR,                           &
             './IOData/'//NULLCHR,                  &
             Debug,VIsdouble)



  ! open(unit=22,file=filename,status='replace')  

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

     xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
     xcell = xedge(:) + dx/2.0;
    
     yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
     ycell = yedge(:) + dy/2.0;
    
     zedge = coord(KAXIS) - bsize(KAXIS)/2.0 + dz*intsz;
     zcell = zedge(:) + dz/2.0; 


     ! Write Block Results into data file:
     write(index_lb,"(I6.6)") lb

     i = TecZne(                                                       &
                'ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR,  &
                 NXB+1,NYB+1,NZB+1,                                    &
                 'BLOCK'//NULLCHR,                                     &
                 CHAR(0))

            
     ! Write x:
     do k=1,NZB+1
        do j=1,NYB+1
           do i=1,NXB+1
              arraylb(i,j,k) = sngl(xedge(i))
           enddo
        enddo
     enddo
     i = TecDat(ijk,arraylb,0)


     ! Write y:
     do k=1,NZB+1
        do j=1,NYB+1
           do i=1,NXB+1
              arraylb(i,j,k) = sngl(yedge(j))
           enddo
        enddo
     enddo
     i = TecDat(ijk,arraylb,0)


     ! Write z:
     do k=1,NZB+1
        do j=1,NYB+1
           do i=1,NXB+1
              arraylb(i,j,k) = sngl(zedge(k))
           enddo
        enddo
     enddo
     i = TecDat(ijk,arraylb,0)

  enddo

  i = TecEnd()

  if (mype .eq. 0) then
  write(*,*) ''
  write(filename,'("./IOData/Grid.",i4.4,".**.plt")') &
        count
  write(*,*) '*** Wrote plotfile to ',filename,' ****'
  endif

  return

66    format(i4.4,g23.15,g23.15,i8.1,i5.1,g23.15)

End subroutine ins_ioWriteTecplot_grid
