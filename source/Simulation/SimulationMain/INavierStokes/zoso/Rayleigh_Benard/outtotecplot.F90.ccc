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
        Grid_releaseBlkPtr, Grid_getBlkIndexLimits, &
        Grid_getBlkBoundBox, Grid_getBlkCenterCoords


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
      
 
  ! Local Variables
  integer :: numblocks,var,i,j,k,lb,nxc,nyc,nzc
  character(29), save :: filename
  character(6) :: index_lb,index_mype

  real xedge(NXB+1),xcell(NXB+1)
  real yedge(NYB+1),ycell(NYB+1)
#ifdef NDIM==3
  real zedge(NZB+1),zcell(NZB+1)
  real intsx(NXB+1),intsy(NYB+1),intsz(NZB+1)
#else
  real intsx(NXB+1),intsy(NYB+1)
#endif

#ifdef NDIM==3
  real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData,facezData
#else
  real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData
#endif

#ifdef NDIM==3
  real facevarxx(NXB+2*NGUARD+1,NYB+2*NGUARD,NZB+2*NGUARD), &
       facevaryy(NXB+2*NGUARD,NYB+2*NGUARD+1,NZB+2*NGUARD), &
       facevarzz(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD+1)
#else
  real facevarxx(NXB+2*NGUARD+1,NYB+2*NGUARD), &
       facevaryy(NXB+2*NGUARD,NYB+2*NGUARD+1), &
#endif

#ifdef NDIM==3
  real, dimension(NXB+1,NYB+1,NZB+1) :: tpu,tpv,tpw,tpp,tpt, &
           tpdudxcorn, tpdudycorn, tpdudzcorn, &
           tpdvdxcorn, tpdvdycorn, tpdvdzcorn, &
           tpdwdxcorn, tpdwdycorn, tpdwdzcorn, &
           vortx, vorty, vortz, divpp, nxp, nyp, nzp
  real*4 arraylb(NXB+1,NYB+1,NZB+1)
  real, dimension(NXB+2*NGUARD,NYB+2*NGUARD, NZB+2*NGUARD) :: tpdudxc, &
        tpdudyc,tpdudzc,tpdvdxc,tpdvdyc,tpdvdzc,tpdwdxc,&
        tpdwdyc,tpdwdzc
#else
  real, dimension(NXB+1,NYB+1,1) :: tpu,tpv,tpp,tpt, &
           tpdudxcorn, tpdudycorn, &
           tpdvdxcorn, tpdvdycorn, &
           vortz, divpp 
  real*4 arraylb(NXB+1,NYB+1,1)
  real, dimension(NXB+2*NGUARD,NYB+2*NGUARD) :: tpdudxc, &
        tpdudyc,tpdvdxc,tpdvdyc
#endif

  integer blockID

  real del(MDIM),dx,dy,dz
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
#ifdef NDIM==3
  ijk       = (NXB+1)*(NYB+1)*(NZB+1)
#else
  ijk       = (NXB+1)*(NYB+1)
#endif
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
#ifdef NDIM==3
  nzc = NZB + NGUARD + 1
#endif

  ! write solution data to data.XXXX.XX
  write(filename,'("./IOData/data.",i4.4,".",i6.6,".plt")') &
        count, mype


#ifdef NDIM==3
  i = TecIni('AMR3D'//NULLCHR,                                      &
             'x y z u v w p t wx wy wz div nx ny nz'//NULLCHR,      &
             filename//NULLCHR,                                     &
             './IOData/'//NULLCHR,                                  &
             Debug,VIsdouble)
#else
  i = TecIni('AMR3D'//NULLCHR,                                      &
             'x y u v p t wz div'//NULLCHR,      &
             filename//NULLCHR,                                     &
             './IOData/'//NULLCHR,                                  &
             Debug,VIsdouble)
#endif

  !open(unit=22,file=filename,status='replace')  

  intsx    = (/ (real(i), i=0,NXB) /)
  intsy    = (/ (real(i), i=0,NYB) /)
#ifdef NDIM==3
  intsz    = (/ (real(i), i=0,NZB) /)
#endif

  call int2char(mype,index_mype)


  do lb = 1,blockcount


     blockID =  blockList(lb)      


     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)
     dx = del(IAXIS)
     dy = del(JAXIS)
#ifdef NDIM==3
     dz = del(KAXIS)
#endif

     ! Get Coord and Bsize for the block:
     ! Bounding box:
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)

     call Grid_getBlkCenterCoords(blockId,coord)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#ifdef NDIM==3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

     tpu = 0.
     tpv = 0.
#ifdef NDIM==3
     tpw = 0.
#endif
     tpp = 0.

     xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
     xcell = xedge(:) + dx/2.0;
    
     yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
     ycell = yedge(:) + dy/2.0;

#ifdef NDIM==3     
     zedge = coord(KAXIS) - bsize(KAXIS)/2.0 + dz*intsz;
     zcell = zedge(:) + dz/2.0;
#endif

     facevarxx = facexData(VELC_FACE_VAR,:,:,:)
     facevaryy = faceyData(VELC_FACE_VAR,:,:,:)
#ifdef NDIM==3
     facevarzz = facezData(VELC_FACE_VAR,:,:,:) 
#endif

#ifdef NDIM==3
     ! U velocity: u(nxb+1,nyb+1,nzb+1)
     ! --------------------------
     tpu = 0.25*(facevarxx(NGUARD+1:nxc,NGUARD:nyc-1,NGUARD:nzc-1) + &
                 facevarxx(NGUARD+1:nxc,NGUARD:nyc-1,NGUARD+1:nzc) + &
                 facevarxx(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc) + &
                 facevarxx(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD:nzc-1));


     ! V velocity: v(nxb+1,nyb+1,nzb+1)
     ! --------------------------                           
     tpv = 0.25*(facevaryy(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD:nzc-1) + &
                 facevaryy(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD:nzc-1) + &
                 facevaryy(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc) + &
                 facevaryy(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD+1:nzc));                              


     ! W velocity: v(nxb+1,nyb+1,nzb+1)
     ! --------------------------                           
     tpw = 0.25*(facevarzz(NGUARD:nxc-1,NGUARD:nyc-1,NGUARD+1:nzc) + &
                 facevarzz(NGUARD+1:nxc,NGUARD:nyc-1,NGUARD+1:nzc) + &
                 facevarzz(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc) + &
                 facevarzz(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD+1:nzc));                              
     

     ! P pressure: p(nxb+1,nyb+1,nzb+1)
     ! -------------------------------
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc, &
                             solnData(PRES_VAR,:,:,:),tpp)

     ! T temperature: t(nxb+1, nyb+1, nzb+1) 
     ! -------------------------------
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc, &
                             solnData(TEMP_VAR,:,:,:),tpt)

     ! Divergence: 
     ! ----------
     solnData(DUST_VAR,:,:,:) = 0.
     solnData(DUST_VAR,NGUARD:nxc,NGUARD:nyc,NGUARD:nzc) =      &
             (facevarxx(NGUARD+1:nxc+1,NGUARD:nyc,NGUARD:nzc) - &
              facevarxx(NGUARD:nxc,NGUARD:nyc,NGUARD:nzc))/dx + &
             (facevaryy(NGUARD:nxc,NGUARD+1:nyc+1,NGUARD:nzc) - &
              facevaryy(NGUARD:nxc,NGUARD:nyc,NGUARD:nzc))/dy + &
             (facevarzz(NGUARD:nxc,NGUARD:nyc,NGUARD+1:nzc+1) - &
              facevarzz(NGUARD:nxc,NGUARD:nyc,NGUARD:nzc))/dz

     solnData(DUST_VAR,1:NGUARD,:,:)=0.
     solnData(DUST_VAR,nxc:2*NGUARD+NXB,:,:)=0.


     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc, &
                             solnData(DUST_VAR,:,:,:),divpp)

     ! Velocity derivatives:
     ! -------- -----------            
     ! Extrapolation of center derivatives to corners, the values
     ! of derivatives in guardcells next to edges are obtained 
     ! from real velocities and linearly extrapolated velocities
     ! to edge points.

     ! U derivatives:
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             tpdudxc,tpdudxcorn)
            
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             tpdudyc,tpdudycorn)

     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             tpdudzc,tpdudzcorn)

     ! V derivatives:
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             tpdvdxc,tpdvdxcorn)
            
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             tpdvdyc,tpdvdycorn)

     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             tpdvdzc,tpdvdzcorn)


     ! W derivatives:
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             tpdwdxc,tpdwdxcorn)
            
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             tpdwdyc,tpdwdycorn)

     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             tpdwdzc,tpdwdzcorn)
         
     ! VORTICITY:
     ! ---------
     ! Corner values of vorticity:
     vortx = tpdwdycorn - tpdvdzcorn
     vorty = tpdudzcorn - tpdwdxcorn
     vortz = tpdvdxcorn - tpdudycorn

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

      ! Write u:
      arraylb(:,:,:) = sngl(tpu)
      i = TecDat(ijk,arraylb,0)

      ! Write v:
      arraylb(:,:,:) = sngl(tpv)
      i = TecDat(ijk,arraylb,0)

      ! Write w:
      arraylb(:,:,:) = sngl(tpw)
      i = TecDat(ijk,arraylb,0)

      ! Write p:
      arraylb(:,:,:) = sngl(tpp)
      i = TecDat(ijk,arraylb,0)

      ! Write t:
      arraylb(:,:,:) = sngl(tpt)
      i = TecDat(ijk,arraylb,0)

      ! Write omgX:
      arraylb(:,:,:) = sngl(vortx)
      i = TecDat(ijk,arraylb,0)

      ! Write omgY:
      arraylb(:,:,:) = sngl(vorty)
      i = TecDat(ijk,arraylb,0)

      ! Write omgZ:
      arraylb(:,:,:) = sngl(vortz)
      i = TecDat(ijk,arraylb,0)

      ! Write Div:
      arraylb(:,:,:) = sngl(divpp)
      i = TecDat(ijk,arraylb,0)

      arraylb(:,:,:) = sngl(nxp)
      i = TecDat(ijk,arraylb,0)

      arraylb(:,:,:) = sngl(nyp)
      i = TecDat(ijk,arraylb,0)

      arraylb(:,:,:) = sngl(nzp)
      i = TecDat(ijk,arraylb,0)

   enddo

#else

     ! U velocity: u(nxb+1,nyb+1)
     ! --------------------------
     tpu = 0.5*(facevarxx(NGUARD+1:nxc,NGUARD:nyc-1) + &
                facevarxx(NGUARD+1:nxc,NGUARD+1:nyc) )


     ! V velocity: v(nxb+1,nyb+1)
     ! --------------------------                           
     tpv = 0.5*(facevaryy(NGUARD:nxc-1,NGUARD+1:nyc) + &
                facevaryy(NGUARD+1:nxc,NGUARD+1:nyc) )


     ! P pressure: p(nxb+1,nyb+1)
     ! -------------------------------
     call centervals2corners(NGUARD,NXB,NYB,nxc,nyc, &
                             solnData(PRES_VAR,:,:,1),tpp)

     ! T temperature: t(nxb+1, nyb+1) 
     ! -------------------------------
     call centervals2corners(NGUARD,NXB,NYB,nxc,nyc, &
                             solnData(TEMP_VAR,:,:,1),tpt)

     ! Divergence: 
     ! ----------
     solnData(DUST_VAR,:,:,:) = 0.
     solnData(DUST_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,1) =      &
             (facevarxx(NGUARD+2:nxc,NGUARD+1:nyc-1) -         &
              facevarxx(NGUARD+1:nxc-1,NGUARD+1:nyc-1))/dx +   &
             (facevaryy(NGUARD+1:nxc-1,NGUARD+2:nyc) -         &
              facevaryy(NGUARD+1:nxc-1,NGUARD+1:nyc-1))/dy

     call centervals2corners(NGUARD,NXB,NYB,nxc,nyc, &
                             solnData(DUST_VAR,:,:,1),divpp)

     ! Velocity derivatives:
     ! -------- -----------            
     ! Extrapolation of center derivatives to corners, the values
     ! of derivatives in guardcells next to edges are obtained 
     ! from real velocities and linearly extrapolated velocities
     ! to edge points.
      tpdudycorn(1:NXB+1,1:NYB+1)=(facevarxx(NGUARD+1:nxc,NGUARD+1:nyc)-  &
                                   facevarxx(NGUARD+1:nxc,NGUARD:nyc-1))/dy

      tpdvdxcorn(1:NXB+1,1:NYB+1)=(facevaryy(NGUARD+1:nxc,NGUARD+1:nyc)-  &
                                   facevaryy(NGUARD:nxc-1,NGUARD+1:nyc))/dx
         
     ! VORTICITY:
     ! ---------
     ! Corner values of vorticity:
     vortz = tpdvdxcorn - tpdudycorn

     ! Write Block Results into data file:
     call int2char(lb,index_lb)

     i = TecZne(                                                      &
                'ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR, &
                 NXB+1,NYB+1,                                         &
                'BLOCK'//NULLCHR,                                     &
                 CHAR(0))

            
      ! Write x:
      do k=1,NZB+1
         do j=1,NYB+1
               arraylb(i,j,1) = sngl(xedge(i))
         enddo
      enddo
      i = TecDat(ijk,arraylb,0)


      ! Write y:
      do k=1,NZB+1
         do j=1,NYB+1
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

      ! Write t:
      arraylb(:,:,1) = sngl(tpt)
      i = TecDat(ijk,arraylb,0)

      ! Write omgZ:
      arraylb(:,:,1) = sngl(vortz)
      i = TecDat(ijk,arraylb,0)

      ! Write Div:
      arraylb(:,:,1) = sngl(divpp)
      i = TecDat(ijk,arraylb,0)

   enddo
#endif

   i = TecEnd()


   if (mype .eq. 0) then
     write(*,*) ''
     write(filename,'("./IOData/data.",i4.4,".**.plt")') &
           count
     write(*,*) '*** Wrote plotfile to ',filename,' ****'
   endif

  End subroutine outtotecplot

  
#ifdef NDIM==3


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

#else

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
