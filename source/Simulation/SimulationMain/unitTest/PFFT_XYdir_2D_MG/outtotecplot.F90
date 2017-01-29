! Subroutine outtotecplot
!
! Subroutine to write out to Tecplot data in binary form.
!
! Written: Marcos Vanella, December 2006.
! ---------------------------------------------------------------------------

subroutine outtotecplot(mype,time,dt,istep,count, &
     timer,blockList,blockCount,firstfileflag)


  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkPtr, &
    Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_getBlkBoundBox, &
    Grid_getBlkCenterCoords


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
  integer :: numblocks,var,i,j,k,lb,nxc,nyc,nzc
  character(25) :: filename
  character(6) :: index_lb,index_mype

  real xedge(NXB+1),xcell(NXB+1)
  real yedge(NYB+1),ycell(NYB+1)
  real intsx(NXB+1), intsy(NYB+1)

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  real facevarxx(NXB+2*NGUARD+1,NYB+2*NGUARD),&
       facevaryy(NXB+2*NGUARD,NYB+2*NGUARD+1),&
       facevarzz(NXB+2*NGUARD,NYB+2*NGUARD)

  real, dimension(NXB+1,NYB+1) :: tpu,tpv,tpw,tpp,&
       tpdudxcorn, tpdudycorn, tpdudzcorn,& 
       tpdvdxcorn, tpdvdycorn, tpdvdzcorn,&
       tpdwdxcorn, tpdwdycorn, tpdwdzcorn,&
       vortx,vorty,vortz,omg,             &
       Sxy,Syz,Sxz,Oxy,Oyz,Oxz,Qcr,divpp,TVtpp


  real*4 arraylb(NXB+1,NYB+1)

  real, dimension(NXB+2*NGUARD,NYB+2*NGUARD) :: tpdudxc,&
       tpdudyc,tpdudzc,tpdvdxc,tpdvdyc,tpdvdzc,tpdwdxc,&
       tpdwdyc,tpdwdzc


  integer blockID

  real del(MDIM),dx,dy
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
  ijk       = (NXB+1)*(NYB+1)
  !-----------------------------------------------------------------------


  ! -- filetime.XX --

!!$  write(filename, '("IOData/data_time.", i2.2)') mype
!!$
!!$  ! create/clear filetime.XX if time = 0
!!$  if(firstfileflag .eq. 0) then
!!$     open(unit=33, file=filename, status='replace')
!!$
!!$#ifdef FLASH_GRID_UG
!!$     write(33,*) 'NXB, NYB, NZB'
!!$     write(33,'(3i4.1)') NXB, NYB, NZB    
!!$#else
!!$     write(33,*) 'NXB, NYB, NZB, interp. order (prolong, restrict)'
!!$     write(33,'(3i4.1)') NXB, NYB, NZB         
!!$#endif
!!$
!!$     write(33,'(a23,a43,a49,a12)') 'file number, time, dt, ',&
!!$          'step number, ',&
!!$          'total number of blocks, number of blocks output, ',&
!!$          'elapsed time'
!!$
!!$     close(33)
!!$  endif
!!$
!!$  ! write timestep data to filetime.XX on each processor
!!$
!!$  open(unit=33, file=filename, status='old', position='append')
!!$  write(33,66) count, time, dt, istep,blockcount,timer
!!$  close(33)
!!$
!!$
!!$  ! -- data.XXXX.XX --
!!$
!!$  nxc = NXB + NGUARD + 1
!!$  nyc = NYB + NGUARD + 1
!!$
!!$  ! write solution data to data.XXXX.XX
!!$
!!$  write(filename,'("./IOData/data.",i4.4,".",i2.2,".plt")') &
!!$       count, mype
!!$
!!$
!!$  i = TecIni('AMR3D'//NULLCHR,                            &
!!$       'x y src anl phi err'//NULLCHR,  &
!!$       filename//NULLCHR,                           &
!!$       './IOData/'//NULLCHR,                  &
!!$       Debug,VIsdouble)
!!$
!!$
!!$
!!$  ! open(unit=22,file=filename,status='replace')  
!!$
!!$  intsx    = (/ (real(i), i=0,NXB) /)
!!$  intsy    = (/ (real(i), i=0,NYB) /)
!!$
!!$  call int2char(mype,index_mype)
!!$
!!$
!!$  do lb = 1,blockcount
!!$
!!$     blockID =  blockList(lb)      
!!$
!!$
!!$     ! Get blocks dx, dy ,dz:
!!$     call Grid_getDeltas(blockID,del)
!!$     dx = del(IAXIS)
!!$     dy = del(JAXIS)
!!$
!!$     ! Get Coord and Bsize for the block:
!!$     ! Bounding box:
!!$     call Grid_getBlkBoundBox(blockId,boundBox)
!!$     bsize(1:NDIM) = boundBox(2,1:NDIM) - boundBox(1,1:NDIM)
!!$
!!$     call Grid_getBlkCenterCoords(blockId,coord)
!!$
!!$     ! Point to blocks center and face vars:
!!$     call Grid_getBlkPtr(blockID,solnData,CENTER)
!!$
!!$     tpu = 0.
!!$     tpv = 0.
!!$     tpp = 0.
!!$     omg = 0.
!!$
!!$     xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
!!$     xcell = xedge(:) + dx/2.0;
!!$
!!$     yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
!!$     ycell = yedge(:) + dy/2.0;
!!$
!!$
!!$     ! source: p(NXB+1,NYB+1,NZB+1)
!!$     ! -------------------------------
!!$     call centervals2corners(NGUARD,NXB,NYB,nxc,nyc, &
!!$          solnData(VSRC_VAR,:,:,1),tpu)
!!$
!!$     ! anl: p(NXB+1,NYB+1,NZB+1)
!!$     ! -------------------------------
!!$     call centervals2corners(NGUARD,NXB,NYB,nxc,nyc, &
!!$          solnData(VANL_VAR,:,:,1),tpp)
!!$
!!$
!!$     ! phi : phi(NXB+1,NYB+1,NZB+1)
!!$     ! -------------------------------
!!$     call centervals2corners(NGUARD,NXB,NYB,nxc,nyc,&
!!$          solnData(VPHI_VAR,:,:,1),TVtpp)
!!$
!!$     ! Error: 
!!$     ! ----------
!!$     call centervals2corners(NGUARD,NXB,NYB,nxc,nyc,&
!!$          solnData(VERR_VAR,:,:,1),divpp)
!!$
!!$
!!$     ! Write Block Results into data file:
!!$     call int2char(lb,index_lb)
!!$
!!$     i = TecZne(                                                       &
!!$          'ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR,  &
!!$          NXB+1,NYB+1,1,                                    &
!!$          'BLOCK'//NULLCHR,                                     &
!!$          CHAR(0))
!!$
!!$
!!$     ! Write x:
!!$     do k=1,1
!!$        do j=1,NYB+1
!!$           do i=1,NXB+1
!!$              arraylb(i,j) = sngl(xedge(i))
!!$           enddo
!!$        enddo
!!$     enddo
!!$     i = TecDat(ijk,arraylb,0)
!!$
!!$
!!$     ! Write y:
!!$     do k=1,1
!!$        do j=1,NYB+1
!!$           do i=1,NXB+1
!!$              arraylb(i,j) = sngl(yedge(j))
!!$           enddo
!!$        enddo
!!$     enddo
!!$     i = TecDat(ijk,arraylb,0)
!!$
!!$
!!$     ! Write source:
!!$     arraylb = sngl(tpu)
!!$     i = TecDat(ijk,arraylb,0)
!!$
!!$     ! Write anl:
!!$     arraylb = sngl(tpp)
!!$     i = TecDat(ijk,arraylb,0)
!!$
!!$     ! Write TV:
!!$     arraylb = sngl(TVtpp)
!!$     i = TecDat(ijk,arraylb,0)
!!$
!!$     ! Write Div:
!!$     arraylb = sngl(divpp)
!!$     i = TecDat(ijk,arraylb,0)
!!$
!!$
!!$  enddo
!!$
!!$  i = TecEnd()
!!$
!!$  return
!!$
!!$66 format(i4.4,g23.15,g23.15,i8.1,i5.1,g23.15)

End subroutine outtotecplot


!!$! Subroutine centervals2corners:
!!$! Subroutine to obtain corver values of a variable given the center 
!!$! values of it in a 2D structured block, suppossing guardcells already
!!$! filled.
!!$!
!!$! ----------------------------------------------------------------------
!!$
!!$       subroutine centervals2corners(ng,nxb,nyb,nxc,nyc,unk1,tpp)
!!$
!!$       implicit none
!!$
!!$       integer ng,nxb,nyb,nxc,nyc
!!$       integer nx1,ny1,nx2,ny2
!!$       real*8, intent(in) :: unk1(nxb+2*ng,nyb+2*ng)
!!$       real*8, intent(out) :: tpp(nxb+1,nyb+1)
!!$
!!$
!!$       tpp = .25*(unk1(ng:nxc-1,ng:nyc-1) + &
!!$                  unk1(ng:nxc-1,ng+1:nyc) + &
!!$                  unk1(ng+1:nxc,ng:nyc-1) + &
!!$                  unk1(ng+1:nxc,ng+1:nyc))
!!$
!!$       ! Corners:
!!$       ! Edge: x = 1, y = 1:
!!$       tpp(1,1) = .5*(unk1(ng,ng+1) + &
!!$                      unk1(ng+1,ng))
!!$
!!$       ! Edge: x = 1, y = end:
!!$       tpp(1,nyb+1) = .5*(unk1(ng,nyc-1) + &
!!$                          unk1(ng+1,nyc))
!!$
!!$       ! Edge: x = end, y = 1:
!!$       tpp(nxb+1,1) = .5*(unk1(nxc-1,ng) + &
!!$                          unk1(nxc,ng+1))
!!$
!!$       ! Edge: x = end, y = end:
!!$       tpp(nxb+1,nyb+1) = .5*(unk1(nxc-1,nyc) + &
!!$                              unk1(nxc,nyc-1))
!!$
!!$
!!$     End subroutine centervals2corners
!!$
!!$! Subroutine int2char
!!$! Subroutine that converts an integer of at most 6 figures
!!$! into a character stored in string
!!$!
!!$! Written by Marcos Vanella in June 2006
!!$! ---------------------------------------------------------------
!!$
!!$Subroutine int2char(i,strng)
!!$
!!$  implicit none
!!$
!!$  integer i
!!$  character (6) strng
!!$
!!$  integer k, val, valaux
!!$  real*8 val2
!!$
!!$  valaux=0 
!!$  strng = '000000'
!!$
!!$
!!$  do k = 6,1,-1
!!$
!!$     val2 = (i-valaux) / (10**(k-1))
!!$     val = floor(val2) 
!!$
!!$     valaux = valaux + val*(10**(k-1))
!!$
!!$     !        write(*,*) 7-k,val,valaux
!!$
!!$     if (val .GE. 1) then
!!$
!!$        select case (val)
!!$
!!$        case (1)
!!$
!!$           strng(7-k:7-k) = "1"
!!$
!!$        case (2)
!!$
!!$           strng(7-k:7-k) = '2'
!!$
!!$
!!$        case (3)
!!$
!!$           strng(7-k:7-k) = '3'
!!$
!!$        case (4)
!!$
!!$           strng(7-k:7-k) = '4'
!!$
!!$        case (5)
!!$
!!$           strng(7-k:7-k) = '5'
!!$
!!$        case (6)
!!$
!!$           strng(7-k:7-k) = '6'
!!$
!!$        case (7)
!!$
!!$           strng(7-k:7-k) = '7'
!!$
!!$        case (8)
!!$
!!$           strng(7-k:7-k) = '8'
!!$
!!$        case (9)
!!$
!!$           strng(7-k:7-k) = '9'
!!$
!!$        end select
!!$
!!$     endif
!!$
!!$  enddo
!!$
!!$End subroutine int2char
!!$
!!$
