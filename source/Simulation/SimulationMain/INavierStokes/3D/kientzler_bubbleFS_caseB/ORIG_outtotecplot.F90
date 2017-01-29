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

!#include "constants.h"
!#include "Flash.h"
  include "Flash_mpi.h"


  integer, intent(in) :: mype,istep,count,firstfileflag
  integer, intent(in) :: blockCount
  integer, intent(in) :: blockList(MAXBLOCKS)
  real, intent(in) :: time,dt,timer
      
 
  ! Local Variables
  integer :: numblocks,var,i,j,k,lb,nxc,nyc,nzc
  character(25) :: filename
  character(6) :: index_lb,index_mype

  real xedge(NXB+1),xcell(NXB+1)
  real yedge(NYB+1),ycell(NYB+1)
  real zedge(NZB+1),zcell(NZB+1)

  real intsx(NXB+1), intsy(NYB+1), intsz(NZB+1)


  real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData,facezData

  real facevarxx(NXB+2*NGUARD+1,NYB+2*NGUARD,NZB+2*NGUARD), &
       facevaryy(NXB+2*NGUARD,NYB+2*NGUARD+1,NZB+2*NGUARD),&
       facevarzz(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD+1)

  real, dimension(NXB+1,NYB+1,NZB+1) :: tpu,tpv,tpw,tpp, &
       tpdudxcorn, tpdudycorn, tpdudzcorn,&
       tpdvdxcorn, tpdvdycorn, tpdvdzcorn,&
       tpdwdxcorn, tpdwdycorn, tpdwdzcorn,&
       vortz,divpp,tpdens,tpdfun,tpvisc

  real*4 arraylb(NXB+1,NYB+1,NZB+1)

  real, dimension(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD) :: tpdudxc,&
       tpdudyc,tpdudzc,tpdvdxc,tpdvdyc,tpdvdzc,tpdwdxc,&
       tpdwdyc,tpdwdzc


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
  ijk       = (NXB+1)*(NYB+1)
!-----------------------------------------------------------------------


! -- filetime.XX --
  write(filename, '("./IOData/data_time.", i2.2)') mype

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
  nzc = NZB + NGUARD + 1

  ! write solution data to data.XXXX.XX
  write(filename,'("./IOData/data.",i4.4,".",i2.2,".plt")') &
        count, mype


  i = TecIni('AMR3D'//NULLCHR,'x y z u v w p dens dfun visc wz div'//NULLCHR,   &
           filename//NULLCHR,'./IOData/'//NULLCHR, &
           Debug,VIsdouble)

  !open(unit=22,file=filename,status='replace')  

  intsx    = (/ (real(i), i=0,NXB) /)
  intsy    = (/ (real(i), i=0,NYB) /)
  intsz    = (/ (real(i), i=0,NZB) /)

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


     tpu = 0.
     tpv = 0.
     tpp = 0.


     xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
     xcell = xedge(:) + dx/2.0;
    
     yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
     ycell = yedge(:) + dy/2.0;

     zedge = coord(KAXIS) - bsize(KAXIS)/2.0 + dz*intsz;
     zcell = zedge(:) + dy/2.0;
    
     facevarxx = facexData(VELC_FACE_VAR,:,:,:)
     facevaryy = faceyData(VELC_FACE_VAR,:,:,:)
     facevarzz = facezData(VELC_FACE_VAR,:,:,:)
 
     ! U velocity: u(nxb+1,nyb+1)
     ! --------------------------
     !tpu = 0.5*(facevarxx(NGUARD+1:nxc,NGUARD:nyc-1)+  &
     !           facevarxx(NGUARD+1:nxc,NGUARD+1:nyc) )
     tpu = 0.0

     ! V velocity: v(nxb+1,nyb+1)
     ! --------------------------                           
     !tpv = 0.5*(facevaryy(NGUARD:nxc-1,NGUARD+1:nyc) + &
     !           facevaryy(NGUARD+1:nxc,NGUARD+1:nyc) )                               
     tpv = 0.0

     ! W velocity: w(nxb+1,nyb+1)
     ! --------------------------                           
     !tpw = 0.5*(facevarzz(NGUARD:nzc-1,NGUARD+1:nzc) + &
     !           facevarzz(NGUARD+1:nzc,NGUARD+1:nzc) )                               
     tpw = 0.0

     ! P pressure: p(nxb+1,nyb+1)
     ! -------------------------------
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                            solnData(PRES_VAR,:,:,1),tpp)

     ! Distance Function: df(nxb+1,nyb+1)
     ! -------------------------------
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                            solnData(DFUN_VAR,:,:,1),tpdfun)

     ! Density: dens(nxb+1,nyb+1)
     ! -------------------------------
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                            solnData(DENS_VAR,:,:,1),tpdens)
     !call centervals2corners(NGUARD,NXB,NYB,nxc,nyc, &
     !                        1./( facexData(RH1F_FACE_VAR,:,:,1) + &
     !                             facexData(RH2F_FACE_VAR,:,:,1) ),tpdens)

     ! Viscosity: visc(nxb+1,nyb+1)
     ! -------------------------------
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                            solnData(VISC_VAR,:,:,1),tpvisc)

     ! Divergence: 
     ! ----------
     solnData(DUST_VAR,NGUARD:nxc,NGUARD:nyc,NGUARD:nzc) =      &
             (facevarxx(NGUARD+1:nxc+1,NGUARD:nyc,NGUARD:nyc) - &
              facevarxx(NGUARD:nxc,NGUARD:nyc,NGUARD:nyc))/dx + &
             (facevaryy(NGUARD:nxc,NGUARD+1:nyc+1,NGUARD:nyc) - &
              facevaryy(NGUARD:nxc,NGUARD:nyc,NGUARD:nyc))/dy + &
             (facevarzz(NGUARD:nxc,NGUARD:nyc,NGUARD+1:nyc+1) - &
              facevarzz(NGUARD:nxc,NGUARD:nyc,NGUARD:nyc))/dz
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             solnData(DUST_VAR,:,:,1),divpp)


            ! Velocity derivatives:
            ! -------- -----------            
!            tpdudxc(ng:nxc,ng:nyc) = (facevarxx(ng+1:nxc+1,ng:nyc) -
!&                                     facevarxx(ng:nxc,ng:nyc))/dx

!            tpdvdyc(ng:nxc,ng:nyc) =  (facevaryy(ng:nxc,ng+1:nyc+1) -
!&                                      facevaryy(ng:nxc,ng:nyc))/dy 

      !tpdudycorn(1:NXB+1,1:NYB+1)=(facevarxx(NGUARD+1:nxc,NGUARD+1:nyc)-  &
      !                             facevarxx(NGUARD+1:nxc,NGUARD:nyc-1))/dy
      tpdudycorn(1:NXB+1,1:NYB+1,1:NZB+1)=0.

      !tpdvdxcorn(1:NXB+1,1:NYB+1)=(facevaryy(NGUARD+1:nxc,NGUARD+1:nyc)-  &
      !                             facevaryy(NGUARD:nxc-1,NGUARD+1:nyc))/dx 
      tpdvdxcorn(1:NXB+1,1:NYB+1,1:NZB+1)=0.
         
      ! VORTICITY:
      ! ---------
      ! Corner values of vorticity:
      vortz = tpdvdxcorn - tpdudycorn


      ! Write Block Results into data file:
      call int2char(lb,index_lb)

      i = TecZne('ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR, &
          NXB+1,NYB+1,1,'BLOCK'//NULLCHR,CHAR(0))

            
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
            arraylb(i,j,k) = sngl(zedge(j))
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

      ! Write dens:
      arraylb(:,:,:) = sngl(tpdens)
      i = TecDat(ijk,arraylb,0)

      ! Write dfun:
      arraylb(:,:,:) = sngl(tpdfun)
      i = TecDat(ijk,arraylb,0)

      ! Write visc:
      arraylb(:,:,:) = sngl(tpvisc)
      i = TecDat(ijk,arraylb,0)


      ! Write omgZ:
      arraylb(:,:,:) = sngl(vortz)
      i = TecDat(ijk,arraylb,0)

      ! Write Div:
      arraylb(:,:,:) = sngl(divpp)
      i = TecDat(ijk,arraylb,0)

   enddo

   i = TecEnd()


  End subroutine outtotecplot

! Subroutine centervals2corners:
! Subroutine to obtain corver values of a variable given the center 
! values of it in a 3D structured block, suppossing guardcells already
! filled.
!
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

! ----------------------------------------------------------------------

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

