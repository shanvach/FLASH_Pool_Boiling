! Subroutine outtotecplot
!
! Subroutine to write out to Tecplot data in binary form.
!
! ---------------------------------------------------------------------------

#include "constants.h"
#include "Flash.h"

  subroutine outtotecplot(mype,time,dt,istep,count,&
           timer,blockList,blockCount,firstfileflag)

      use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkPtr, &
        Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_getCellCoords, &
        Grid_getBlkBoundBox, Grid_getBlkCenterCoords, Grid_getCellMetrics

      use Driver_data, only : dr_globalMe
#ifdef FLASH_GRID_UG
#else
      use physicaldata, ONLY : interp_mask_unk,interp_mask_unk_res
#endif

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

  real xedge(GRID_IHI+1), xcell(GRID_IHI)
  real yedge(GRID_JHI+1), ycell(GRID_JHI)
  real intsx(NXB+1), intsy(NYB+1)


  real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData

  real facevarxx(NXB+2*NGUARD+1,NYB+2*NGUARD), &
       facevaryy(NXB+2*NGUARD,NYB+2*NGUARD+1), &
       forcvarxx(NXB+2*NGUARD+1,NYB+2*NGUARD), &
       forcvaryy(NXB+2*NGUARD,NYB+2*NGUARD+1), &
       facevarxn(NXB+2*NGUARD+1,NYB+2*NGUARD), &
       facevaryn(NXB+2*NGUARD,NYB+2*NGUARD+1)
 


  real, dimension(NXB+1,NYB+1) :: tpu,tpv,tpp, &
           tpfu,tpfv,tpun,tpvn, &
           tpdudxcorn, tpdudycorn, &
           tpdvdxcorn, tpdvdycorn, &
           vortz,divpp,tpcurv,tpdfun, tpt ! tpt declared by Akash

  real*4 arraylb(NXB+1,NYB+1,1)
 
  real, dimension(NXB+2*NGUARD,NYB+2*NGUARD) :: tpdudxc, &
        tpdudyc,tpdvdxc,tpdvdyc

  real, dimension(GRID_IHI_GC,3,blockCount) :: dx
  real, dimension(GRID_JHI_GC,3,blockCount) :: dy

  integer blockID, ix, jy, locx, locy

  real del(MDIM), mval
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
  !write(*,*) 'In outtotecplot before write mype',mype
  !write(filename, '("./IOData/data_time1234.", i2.2)') mype
  !write(*,*) 'filename=',filename
  !write(*, '("./IOData/data_time.", i2.2)') mype  
  !write(filename, '("./IOData/data_time.", i2.2)') mype
  write(filename, '("./IOData/data_time.", i6.6)') mype
  !write(*,*) 'after write mype',mype

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

  !write(*,*) "Write 2D TecPlot Data"
  ! write solution data to data.XXXX.XX
  !write(filename,'("./IOData/data.",i4.4,".",i2.2,".plt")') &
  write(filename,'("./IOData/data.",i4.4,".",i6.6,".plt")') &
        count, mype


  i = TecIni('AMR2D'//NULLCHR,'x y u v p forcu forcv un vn'//NULLCHR,   &
           filename//NULLCHR,'./IOData/'//NULLCHR, &
           Debug,VIsdouble)

  open(unit=22,file=filename,status='replace')  

  intsx    = (/ (real(i), i=0,NXB) /)
  intsy    = (/ (real(i), i=0,NYB) /)

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
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

  ! Get blk cell metrics by direction from Grid Unit 
    
     call Grid_getCellMetrics(IAXIS,blockID,LEFT_EDGE, .true.,dx(:,LEFT_EDGE,lb), GRID_IHI_GC) 
     call Grid_getCellMetrics(IAXIS,blockID,CENTER,    .true.,dx(:,CENTER,lb),    GRID_IHI_GC) 
     call Grid_getCellMetrics(IAXIS,blockID,RIGHT_EDGE,.true.,dx(:,RIGHT_EDGE,lb),GRID_IHI_GC) 

     call Grid_getCellMetrics(JAXIS,blockID,LEFT_EDGE, .true.,dy(:,LEFT_EDGE,lb), GRID_JHI_GC) 
     call Grid_getCellMetrics(JAXIS,blockID,CENTER,    .true.,dy(:,CENTER,lb),    GRID_JHI_GC) 
     call Grid_getCellMetrics(JAXIS,blockID,RIGHT_EDGE,.true.,dy(:,RIGHT_EDGE,lb),GRID_JHI_GC) 

     tpu = 0.
     tpv = 0.
     tpun = 0.
     tpvn = 0. 
     tpfu = 0.
     tpfv = 0. 
     tpp = 0.
     tpt = 0.
     divpp = 0.
     tpcurv = 0.
     tpdfun = 0.

     call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC)
     call Grid_getCellCoords(IAXIS, blockId, CENTER, .false., xcell, blkLimits(HIGH, IAXIS)) 
     call Grid_getCellCoords(JAXIS, blockId, CENTER, .false., ycell, blkLimits(HIGH, JAXIS)) 
     call Grid_getCellCoords(IAXIS, blockId, FACES, .false., xedge, blkLimits(HIGH, IAXIS)+1)   
     call Grid_getCellCoords(JAXIS, blockId, FACES, .false., yedge, blkLimits(HIGH, JAXIS)+1)   
     !call Grid_getCellCoords(IAXIS, blockId, FACES, .false., xedge, blkLimits(HIGH, IAXIS)+1)  
     !do j=1,blkLimitsGC(HIGH,JAXIS)
        !if(1/dy(j,CENTER,lb) .lt. 0.3125) print*,dy(j,CENTER,lb),j
     !enddo
     !stop

     !print*,maxval(xedge),maxval(yedge)
     !stop
 
     facevarxx = facexData(VELC_FACE_VAR,:,:,1)
     facevaryy = faceyData(VELC_FACE_VAR,:,:,1)
     forcvarxx = facexData(FORC_FACE_VAR,:,:,1)
     forcvaryy = faceyData(FORC_FACE_VAR,:,:,1)
     facevarxn = facexData(RHDS_FACE_VAR,:,:,1)
     facevaryn = faceyData(RHDS_FACE_VAR,:,:,1)
 
     ! U velocity: u(nxb+1,nyb+1)
     ! --------------------------
     tpu = 0.5*(facevarxx(NGUARD+1:nxc,NGUARD:nyc-1)+  &
                facevarxx(NGUARD+1:nxc,NGUARD+1:nyc) )


     ! V velocity: v(nxb+1,nyb+1)
     ! --------------------------                           
     tpv = 0.5*(facevaryy(NGUARD:nxc-1,NGUARD+1:nyc) + &
                facevaryy(NGUARD+1:nxc,NGUARD+1:nyc) )                               


     ! U-new velocity: u(nxb+1,nyb+1)
     ! --------------------------
     tpun = 0.5*(facevarxn(NGUARD+1:nxc,NGUARD:nyc-1)+  &
                facevarxn(NGUARD+1:nxc,NGUARD+1:nyc) )


     ! V-new velocity: v(nxb+1,nyb+1)
     ! --------------------------                           
     tpvn = 0.5*(facevaryn(NGUARD:nxc-1,NGUARD+1:nyc) + &
                facevaryn(NGUARD+1:nxc,NGUARD+1:nyc) )                               



     ! U force: u(nxb+1,nyb+1)
     ! --------------------------
     tpfu = 0.5*(forcvarxx(NGUARD+1:nxc,NGUARD:nyc-1)+  &
                forcvarxx(NGUARD+1:nxc,NGUARD+1:nyc) )


     ! V force: v(nxb+1,nyb+1)
     ! --------------------------                           
     tpfv = 0.5*(forcvaryy(NGUARD:nxc-1,NGUARD+1:nyc) + &
                forcvaryy(NGUARD+1:nxc,NGUARD+1:nyc) )                               



     !print*,size(solnData(:,:,:,:),1),size(solnData(:,:,:,:),4) 
     ! P pressure: p(nxb+1,nyb+1)
     ! -------------------------------
     !print*,NXB+1,GRID_IHI+1,blkLimits(HIGH, IAXIS)+1
     call centervals2corners(NGUARD,NXB,NYB,nxc,nyc, &
                            solnData(PRES_VAR,:,:,1),tpp)
 
     ! T temperature: t(nxb+1, nyb+1) - ! Akash
     ! -------------------------------
     !call centervals2corners(NGUARD,NXB,NYB,nxc,nyc, &
     !                       solnData(TEMP_VAR,:,:,1),tpt)

     ! Divergence
     ! --------------------------------
     mval = 1e-19
     !do ix=1,nxc+NGUARD-1
     !   do jy=1,nyc+NGUARD-1
     !       solnData(DUST_VAR,ix,jy,1) =      &
     !               (facevarxx(ix+1,jy  ) - &
     !                facevarxx(ix  ,jy  ))*dx(ix,CENTER,blockID) + &
     !               (facevaryy(ix  ,jy+1) - &
     !                facevaryy(ix  ,jy  ))*dy(jy,CENTER,blockID)
            !if(solnData(DUST_VAR,ix,jy,1).gt.mval) then
            !   mval = solnData(DUST_VAR,ix,jy,1)
            !   locx = ix
            !   locy = jy
            ! endif

      !open(1, file = 'divu2.txt', action = 'write', access = 'append') !status ='replace' 
      !        if(solnData(DUST_VAR,ix,jy,1).gt.1e-6) write(1,*) solnData(DUST_VAR,ix,jy,1),ix,jy,dr_globalMe
      !close(1)


     !   enddo
     !enddo
 
     !print*,"mvalue",mval,locx,locy
     !print*,"maxdiv",maxval(solnData(DUST_VAR,:,:,1))
     !print*,"mindiv",minval(solnData(DUST_VAR,:,:,1))

     !call centervals2corners(NGUARD,NXB,NYB,nxc,nyc, &
     !                        solnData(DUST_VAR,:,:,1),divpp)
     
    ! Curvature
     ! -------------------------------
     !call centervals2corners(NGUARD,NXB,NYB,nxc,nyc, &
     !                       solnData(CURV_VAR,:,:,1),tpcurv)

    
     ! Distance Function
     ! --------------------------------
     !call centervals2corners(NGUARD,NXB,NYB,nxc,nyc, &
     !                       solnData(DFUN_VAR,:,:,1),tpdfun)

      ! Write Block Results into data file:
      call int2char(lb,index_lb)

      i = TecZne('ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR, &
          NXB+1,NYB+1,1,'BLOCK'//NULLCHR,CHAR(0))

            
      ! Write x:
      do j=1,NYB+1
         do i=1,NXB+1
            arraylb(i,j,1) = sngl(xedge(i))
            !print*,sngl(xedge(i))
         enddo
      enddo
      !print*,ijk,size(arraylb)
      !stop
      i = TecDat(ijk,arraylb,0)


      ! Write y:
      do j=1,NYB+1
         do i=1,NXB+1
            arraylb(i,j,1) = sngl(yedge(j))
         enddo
      enddo
      i = TecDat(ijk,arraylb,0)


      ! Write u:
      arraylb(:,:,1) = sngl(tpu)
      i = TecDat(ijk,arraylb,0)

      ! Write v:
      arraylb(:,:,1) = sngl(tpv)
      i = TecDat(ijk,arraylb,0)

      ! Write p:
      arraylb(:,:,1) = sngl(tpp)
      i = TecDat(ijk,arraylb,0)

      ! Write u force:
      arraylb(:,:,1) = sngl(tpfu)
      i = TecDat(ijk,arraylb,0)

      ! Write v force:
      arraylb(:,:,1) = sngl(tpfv)
      i = TecDat(ijk,arraylb,0)

      ! Write un:
      arraylb(:,:,1) = sngl(tpun)
      i = TecDat(ijk,arraylb,0)

      ! Write vn:
      arraylb(:,:,1) = sngl(tpvn)
      i = TecDat(ijk,arraylb,0)



      ! Write dfun:
      !arraylb(:,:,1) = sngl(tpdfun)
      !i = TecDat(ijk,arraylb,0)

      ! Write curv:
      !arraylb(:,:,1) = sngl(tpcurv)
      !i = TecDat(ijk,arraylb,0)

      ! Write div:
      !arraylb(:,:,1) = sngl(divpp)
      !i = TecDat(ijk,arraylb,0)

      ! Write t: ! Akash
      !arraylb(:,:,1) = sngl(tpt)
      !i = TecDat(ijk,arraylb,0)

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
