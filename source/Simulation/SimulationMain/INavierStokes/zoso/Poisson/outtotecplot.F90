! Subroutine outtotecplot
!
! Subroutine to write out to Tecplot data in binary form.
!
! ---------------------------------------------------------------------------

#include "constants.h"
#include "Flash.h"


#if NDIM == 3
  subroutine outtotecplot(mype,time,dt,istep,count,&
           timer,blockList,blockCount,firstfileflag)

      use Grid_interface, ONLY : Grid_getCellMetrics, Grid_getBlkPtr,   &
        Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_getCellCoords, &
        Grid_getBlkBoundBox, Grid_getBlkCenterCoords

  implicit none

  include "Flash_mpi.h"

  integer, intent(in) :: mype,istep,count,firstfileflag
  integer, intent(in) :: blockCount
  integer, intent(in) :: blockList(MAXBLOCKS)
  real, intent(in) :: time,dt,timer
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
 
  ! Local Variables
  integer :: numblocks,var,i,j,k,lb,nxc,nyc,nzc
  character(29), save :: filename
  character(6) :: index_lb,index_mype

  real xedge(NXB+1),xcell(NXB)
  real yedge(NYB+1),ycell(NYB)
  real zedge(NZB+1),zcell(NZB)
  real intsx(NXB+1),intsy(NYB+1),intsz(NZB+1)

  real, pointer, dimension(:,:,:,:) :: solnData

  real, dimension(GRID_IHI_GC,3,blockCount) :: iMetrics
  real, dimension(GRID_JHI_GC,3,blockCount) :: jMetrics
  real, dimension(GRID_KHI_GC,3,blockCount) :: kMetrics

  real, dimension(NXB+1,NYB+1,NZB+1) :: nsrc,nfld,     &
                                        nxp, nyp, nzp

  real*4 arraylb(NXB+1,NYB+1,NZB+1)

  integer blockID

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


! -- filetime.XX --
  
  write(filename, '("./IOData/data_time.", i4.4)') mype

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

66    format(i4.4,g23.15,g23.15,i8.1,i5.1,g23.15)


  ! -- data.XXXX.XX --
  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1
  nzc = NZB + NGUARD + 1

  ! write(*,*) "Write 3D Tecplot Data"
  ! write solution data to data.XXXX.XX
  write(filename,'("./IOData/data.",i4.4,".",i6.6,".plt")') &
        count, mype


  i = TecIni('AMR3D'//NULLCHR,                                      &
             'x y z nsrc nfld'//NULLCHR,      &
             filename//NULLCHR,                                     &
             './IOData/'//NULLCHR,                                  &
             Debug,VIsdouble)

  !open(unit=22,file=filename,status='replace')  

  intsx    = (/ (real(i), i=0,NXB) /)
  intsy    = (/ (real(i), i=0,NYB) /)
  intsz    = (/ (real(i), i=0,NZB) /)

  call int2char(mype,index_mype)


  do lb = 1,blockcount


     blockID =  blockList(lb)      

     ! Get Coord and Bsize for the block:
     ! Bounding box:
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)

     call Grid_getBlkCenterCoords(blockId,coord)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC)
     call Grid_getCellCoords(IAXIS, blockId, CENTER, .false., xcell, blkLimits(HIGH, IAXIS)-1)
     call Grid_getCellCoords(JAXIS, blockId, CENTER, .false., ycell, blkLimits(HIGH, JAXIS)-1)
     call Grid_getCellCoords(KAXIS, blockId, CENTER, .false., zcell, blkLimits(HIGH, KAXIS)-1)
     call Grid_getCellCoords(IAXIS, blockId, FACES, .false., xedge, blkLimits(HIGH, IAXIS))
     call Grid_getCellCoords(JAXIS, blockId, FACES, .false., yedge, blkLimits(HIGH, JAXIS))
     call Grid_getCellCoords(KAXIS, blockId, FACES, .false., zedge, blkLimits(HIGH, KAXIS))

     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc, &
                             solnData(NSRC_VAR,:,:,:),nsrc)

     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc, &
                             solnData(NFLD_VAR,:,:,:),nfld)

     ! Write Block Results into data file:
     call int2char(lb,index_lb)

     i = TecZne(                                                      &
                'ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR, &
                 NXB+1,NYB+1,NZB+1,                                   &
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

      arraylb(:,:,:) = sngl(nsrc)
      i = TecDat(ijk,arraylb,0)

      arraylb(:,:,:) = sngl(nfld)
      i = TecDat(ijk,arraylb,0)

   enddo

   i = TecEnd()

   if (mype .eq. 0) then
     write(*,*) ''
     write(filename,'("./IOData/data.",i4.4,".**.plt")') &
           count
     write(*,*) '*** Wrote plotfile to ',filename,' ****'
   endif

  End subroutine outtotecplot


! Subroutine centervals2corners:
! Subroutine to obtain corver values of a variable given the center 
! values of it in a 3D structured block, suppossing guardcells already
! filled.
!
! Written by Marcos Vanella Decemeber 2006.
! ----------------------------------------------------------------------

  subroutine centervals2corners(ng,nxb,nyb,nzb,nxc,nyc,nzc, &
                                unk1,tpp)

    implicit none

    integer ng,nxb,nyb,nzb,nxc,nyc,nzc
    integer nx1,ny1,nz1,nx2,ny2,nz2
    real*8, intent(in) :: unk1(nxb+2*ng,nyb+2*ng,nzb+2*ng)
    real*8, intent(out) :: tpp(nxb+1,nyb+1,nzb+1)

      
    tpp = .5*.25*(unk1(ng:nxc-1,ng:nyc-1,ng:nzc-1) + &
                  unk1(ng:nxc-1,ng:nyc-1,ng+1:nzc) + &
                  unk1(ng:nxc-1,ng+1:nyc,ng:nzc-1) + &
                  unk1(ng+1:nxc,ng:nyc-1,ng:nzc-1) + &
                  unk1(ng+1:nxc,ng+1:nyc,ng:nzc-1) + &
                  unk1(ng:nxc-1,ng+1:nyc,ng+1:nzc) + &
                  unk1(ng+1:nxc,ng:nyc-1,ng+1:nzc) + &
                  unk1(ng+1:nxc,ng+1:nyc,ng+1:nzc));

    ! Z edges:
    ! Edge: x = 1, y = 1:
    tpp(1,1,:) = .25*(unk1(ng,ng+1,ng:nzc-1) + &
                      unk1(ng,ng+1,ng+1:nzc) + &
                      unk1(ng+1,ng,ng:nzc-1) + &
                      unk1(ng+1,ng,ng+1:nzc));           
            
    ! Edge: x = 1, y = end:
    tpp(1,nyb+1,:) = .25*(unk1(ng,nyc-1,ng:nzc-1) + & 
                          unk1(ng,nyc-1,ng+1:nzc) + &
                          unk1(ng+1,nyc,ng:nzc-1) + &
                          unk1(ng+1,nyc,ng+1:nzc));    

    ! Edge: x = end, y = 1:
    tpp(nxb+1,1,:) = .25*(unk1(nxc-1,ng,ng:nzc-1) + &
                          unk1(nxc-1,ng,ng+1:nzc) + &
                          unk1(nxc,ng+1,ng:nzc-1) + &
                          unk1(nxc,ng+1,ng+1:nzc));

    ! Edge: x = end, y = end:
    tpp(nxb+1,nyb+1,:) = .25*(unk1(nxc-1,nyc,ng:nzc-1) + &
                              unk1(nxc-1,nyc,ng+1:nzc) + &
                              unk1(nxc,nyc-1,ng:nzc-1) + &
                              unk1(nxc,nyc-1,ng+1:nzc));
            
    ! Y edges:
    ! Edge: x = 1, z = 1:
    tpp(1,:,1) = .25*(unk1(ng,ng:nyc-1,ng+1) + &
                      unk1(ng,ng+1:nyc,ng+1) + &
                      unk1(ng+1,ng:nyc-1,ng) + &
                      unk1(ng+1,ng+1:nyc,ng));

    ! Edge: x = 1, z = end:
    tpp(1,:,nzb+1) = .25*(unk1(ng,ng:nyc-1,nzc-1) + &
                          unk1(ng,ng+1:nyc,nzc-1) + &
                          unk1(ng+1,ng:nyc-1,nzc) + &
                          unk1(ng+1,ng+1:nyc,nzc)); 

    ! Edge: x = end, z = 1:
    tpp(nxb+1,:,1) = .25*(unk1(nxc-1,ng:nyc-1,ng) + &
                          unk1(nxc-1,ng+1:nyc,ng) + &
                          unk1(nxc,ng:nyc-1,ng+1) + &
                          unk1(nxc,ng+1:nyc,ng+1));

    ! Edge: x = end, z = end:
    tpp(nxb+1,:,nzb+1) = .25*(unk1(nxc-1,ng:nyc-1,nzc) + &
                              unk1(nxc-1,ng+1:nyc,nzc) + &
                              unk1(nxc,ng:nyc-1,nzc-1) + &
                              unk1(nxc,ng+1:nyc,nzc-1));

    ! X edges:
    ! Edge: y = 1, z = 1:
    tpp(:,1,1) = .25*(unk1(ng:nxc-1,ng,ng+1) + &
                      unk1(ng+1:nxc,ng,ng+1) + &
                      unk1(ng:nxc-1,ng+1,ng) + &
                      unk1(ng+1:nxc,ng+1,ng));            

    ! Edge: y = 1, z = end:
    tpp(:,1,nzb+1) = .25*(unk1(ng:nxc-1,ng,nzc-1) + &
                          unk1(ng+1:nxc,ng,nzc-1) + &
                          unk1(ng:nxc-1,ng+1,nzc) + &
                          unk1(ng+1:nxc,ng+1,nzc));   

    ! Edge: y = end, z = 1:
    tpp(:,nyb+1,1) = .25*(unk1(ng:nxc-1,nyc-1,ng) + & 
                          unk1(ng+1:nxc,nyc-1,ng) + &
                          unk1(ng:nxc-1,nyc,ng+1) + &
                          unk1(ng+1:nxc,nyc,ng+1)); 

    ! Edge: y = end, z = end:
    tpp(:,nyb+1,nzb+1) = .25*(unk1(ng:nxc-1,nyc-1,nzc) + &
                              unk1(ng+1:nxc,nyc-1,nzc) + &
                              unk1(ng:nxc-1,nyc,nzc-1) + &
                              unk1(ng+1:nxc,nyc,nzc-1)); 

    ! Corners
    ! Corner x = 1, y = 1, z = 1:
    tpp(1,1,1) = -0.5*unk1(ng+1,ng+1,ng+1) + &
                  0.5*unk1(ng,ng+1,ng+1)   + &
                  0.5*unk1(ng+1,ng,ng+1)   + &
                  0.5*unk1(ng+1,ng+1,ng)   


    ! Corner x = end, y =1, z = 1:
    tpp(nxb+1,1,1) = -0.5*unk1(nxc-1,ng+1,ng+1) + &
                      0.5*unk1(nxc,ng+1,ng+1)   + &
                      0.5*unk1(nxc-1,ng,ng+1)   + &
                      0.5*unk1(nxc-1,ng+1,ng)   


    ! Corner x = end, y = end, z = 1:
    tpp(nxb+1,nyb+1,1) = -0.5*unk1(nxc-1,nyc-1,ng+1) + &
                          0.5*unk1(nxc,nyc-1,ng+1)   + &
                          0.5*unk1(nxc-1,nyc,ng+1)   + &
                          0.5*unk1(nxc-1,nyc-1,ng)

    ! Corner x = 1, y = end, z = 1:
    tpp(1,nyb+1,1) = -0.5*unk1(ng+1,nyc-1,ng+1) + &
                      0.5*unk1(ng,nyc-1,ng+1)   + &
                      0.5*unk1(ng+1,nyc,ng+1)   + &
                      0.5*unk1(ng+1,nyc-1,ng)   

    ! Corner x = 1, y = 1, z = end:
    tpp(1,1,nzb+1) = -0.5*unk1(ng+1,ng+1,nzc-1) + &
                      0.5*unk1(ng,ng+1,nzc-1)   + &
                      0.5*unk1(ng+1,ng,nzc-1)   + &
                      0.5*unk1(ng+1,ng+1,nzc)               

    ! Corner x = end, y = 1, z = end:
    tpp(nxb+1,1,nzb+1) = -0.5*unk1(nxc-1,ng+1,nzc-1) + &
                          0.5*unk1(nxc,ng+1,nzc-1)   + &
                          0.5*unk1(nxc-1,ng,nzc-1)   + &
                          0.5*unk1(nxc-1,ng+1,nzc)    

    ! Corner x = end, y = end, z = end:
    tpp(nxb+1,nyb+1,nzb+1) = -0.5*unk1(nxc-1,nyc-1,nzc-1) + &
                              0.5*unk1(nxc,nyc-1,nzc-1)   + &
                              0.5*unk1(nxc-1,nyc,nzc-1)   + &
                              0.5*unk1(nxc-1,nyc-1,nzc)

    ! Corner x = 1, y = end, z = end:
    tpp(1,nyb+1,nzb+1) = -0.5*unk1(ng+1,nyc-1,nzc-1) + &
                          0.5*unk1(ng,nyc-1,nzc-1)   + &
                          0.5*unk1(ng+1,nyc,nzc-1)   + &
                          0.5*unk1(ng+1,nyc-1,nzc)          

  End subroutine centervals2corners

#endif

! Subroutine outtotecplot
!
! Subroutine to write out to Tecplot data in binary form.
!
! ---------------------------------------------------------------------------

!#include "constants.h"
!#include "Flash.h"

#if NDIM == 2

  subroutine outtotecplot(mype,time,dt,istep,count,&
           timer,blockList,blockCount,firstfileflag)

      use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkPtr, &
        Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_getCellCoords, &
        Grid_getBlkBoundBox, Grid_getBlkCenterCoords


  implicit none

  include "Flash_mpi.h"

  integer, intent(in) :: mype,istep,count,firstfileflag
  integer, intent(in) :: blockCount
  integer, intent(in) :: blockList(MAXBLOCKS)
  real, intent(in) :: time,dt,timer
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC    
 
  ! Local Variables
  integer :: numblocks,var,i,j,k,lb,nxc,nyc,nzc
  !character(25), save :: filename
  character(29), save :: filename
  character(6) :: index_lb,index_mype

  real xedge(NXB+1), xcell(NXB)
  real yedge(NYB+1), ycell(NYB)
  real intsx(NXB+1), intsy(NYB+1)

  real, pointer, dimension(:,:,:,:) :: solnData

  real, dimension(NXB+1,NYB+1) :: nsrc,nfld 

  real*4 arraylb(NXB+1,NYB+1,1)
 
  integer blockID

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
  write(filename, '("./IOData/data_time.", i6.6)') mype

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
66    format(i4.4,g23.15,g23.15,i8.1,i5.1,g23.15)


  ! -- data.XXXX.XX --
  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1

  ! write(*,*) "Write 2D TecPlot Data"
  ! write solution data to data.XXXX.XX
  !write(filename,'("./IOData/data.",i4.4,".",i2.2,".plt")') &
  write(filename,'("./IOData/data.",i4.4,".",i6.6,".plt")') &
        count, mype


  i = TecIni('AMR2D'//NULLCHR,'x y nsrc nfld'//NULLCHR,   &
           filename//NULLCHR,'./IOData/'//NULLCHR, &
           Debug,VIsdouble)

  open(unit=22,file=filename,status='replace')  

  intsx    = (/ (real(i), i=0,NXB) /)
  intsy    = (/ (real(i), i=0,NYB) /)

  call int2char(mype,index_mype)


  do lb = 1,blockcount


     blockID =  blockList(lb)      

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC)
     call Grid_getCellCoords(IAXIS, blockId, CENTER, .false., xcell, blkLimits(HIGH, IAXIS)-1) 
     call Grid_getCellCoords(JAXIS, blockId, CENTER, .false., ycell, blkLimits(HIGH, JAXIS)-1) 
     call Grid_getCellCoords(IAXIS, blockId, FACES, .false., xedge, blkLimits(HIGH, IAXIS)) 
     call Grid_getCellCoords(JAXIS, blockId, FACES, .false., yedge, blkLimits(HIGH, JAXIS)) 

     call centervals2corners(NGUARD,NXB,NYB,nxc,nyc, &
                            solnData(NSRC_VAR,:,:,1),nsrc)
     call centervals2corners(NGUARD,NXB,NYB,nxc,nyc, &
                            solnData(NFLD_VAR,:,:,1),nfld)


      ! Write Block Results into data file:
      call int2char(lb,index_lb)

      i = TecZne('ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR, &
          NXB+1,NYB+1,1,'BLOCK'//NULLCHR,CHAR(0))

            
      ! Write x:
      do j=1,NYB+1
         do i=1,NXB+1
            arraylb(i,j,1) = sngl(xedge(i))
         enddo
      enddo
      i = TecDat(ijk,arraylb,0)


      ! Write y:
      do j=1,NYB+1
         do i=1,NXB+1
            arraylb(i,j,1) = sngl(yedge(j))
         enddo
      enddo
      i = TecDat(ijk,arraylb,0)

      arraylb(:,:,1) = sngl(nsrc)
      i = TecDat(ijk,arraylb,0)

      arraylb(:,:,1) = sngl(nfld)
      i = TecDat(ijk,arraylb,0)


   enddo

   i = TecEnd()


  if (mype .eq. 0) then
  write(*,*) ''
  write(filename,'("./IOData/data.",i4.4,".**.plt")') &
        count
  write(*,*) '*** Wrote plotfile to ',filename,' ****'
  endif

  End subroutine outtotecplot


! Subroutine centervals2corners:
! Subroutine to obtain corver values of a variable given the center 
! values of it in a 2D structured block, suppossing guardcells already
! filled.
!
! ----------------------------------------------------------------------

       subroutine centervals2corners(ng,nxb,nyb,nxc,nyc,unk1,tpp)

       implicit none

       integer ng,nxb,nyb,nxc,nyc
       integer nx1,ny1,nx2,ny2
       real*8, intent(in) :: unk1(nxb+2*ng,nyb+2*ng)
       real*8, intent(out) :: tpp(nxb+1,nyb+1)


       tpp = .25*(unk1(ng:nxc-1,ng:nyc-1) + &
                  unk1(ng:nxc-1,ng+1:nyc) + &
                  unk1(ng+1:nxc,ng:nyc-1) + &
                  unk1(ng+1:nxc,ng+1:nyc))

       ! Corners:
       ! Edge: x = 1, y = 1:
       tpp(1,1) = .5*(unk1(ng,ng+1) + &
                      unk1(ng+1,ng))

       ! Edge: x = 1, y = end:
       tpp(1,nyb+1) = .5*(unk1(ng,nyc-1) + &
                          unk1(ng+1,nyc))

       ! Edge: x = end, y = 1:
       tpp(nxb+1,1) = .5*(unk1(nxc-1,ng) + &
                          unk1(nxc,ng+1))

       ! Edge: x = end, y = end:
       tpp(nxb+1,nyb+1) = .5*(unk1(nxc-1,nyc) + &
                              unk1(nxc,nyc-1))


     End subroutine centervals2corners



#endif

! Subroutine int2char
! Subroutine that converts an integer of at most 6 figures
! into a character stored in string
!
! Written by Marcos Vanella in June 2006
! ---------------------------------------------------------------

      Subroutine int2char(i,strng)

      integer i
      character (6) strng
      
      integer k, val, valaux
      real*8 val2

      valaux=0 
      strng = '000000'


      do k = 6,1,-1

         val2 = (i-valaux) / (10**(k-1))
         val = floor(val2) 

         valaux = valaux + val*(10**(k-1))

!        write(*,*) 7-k,val,valaux

         if (val .GE. 1) then

            select case (val)

            case (1)
               
               strng(7-k:7-k) = "1"
               
            case (2)
               
               strng(7-k:7-k) = '2'
               
               
            case (3)
               
               strng(7-k:7-k) = '3'
               
            case (4)
               
               strng(7-k:7-k) = '4'
               
            case (5)
               
               strng(7-k:7-k) = '5'
               
            case (6)
               
               strng(7-k:7-k) = '6'
               
            case (7)
               
               strng(7-k:7-k) = '7'
               
            case (8)
               
               strng(7-k:7-k) = '8'
               
            case (9)
               
               strng(7-k:7-k) = '9'
               
            end select
            
         endif

      enddo

      End subroutine int2char
