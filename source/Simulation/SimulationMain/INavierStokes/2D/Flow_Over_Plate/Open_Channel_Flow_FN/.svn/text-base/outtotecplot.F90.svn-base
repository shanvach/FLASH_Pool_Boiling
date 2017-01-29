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

  use Driver_data,      ONLY : dr_nstep


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
  !character(25) :: filename
  character(29) :: filename
  character(6) :: index_lb,index_mype

  real xedge(NXB+1),xcell(NXB+1)
  real yedge(NYB+1),ycell(NYB+1)
  real intsx(NXB+1), intsy(NYB+1)


  real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData

  real facevarxx(NXB+2*NGUARD+1,NYB+2*NGUARD), &
       facevaryy(NXB+2*NGUARD,NYB+2*NGUARD+1)

  real facevarr1(NXB+2*NGUARD+1,NYB+2*NGUARD), &
       facevarr2(NXB+2*NGUARD,NYB+2*NGUARD+1)
  real facevarr3(NXB+2*NGUARD+1,NYB+2*NGUARD), &
       facevarr4(NXB+2*NGUARD,NYB+2*NGUARD+1)

  real, dimension(NXB+1,NYB+1) :: tpu,tpv,tpp, &
           tpdudxcorn, tpdudycorn, &
           tpdvdxcorn, tpdvdycorn, &
           vortz,divpp,tpdens,tpdensy,tpdfun,tpvisc,tpcurv,qcriterion

  real*4 arraylb(NXB+1,NYB+1,1)
 
  real, dimension(NXB+2*NGUARD,NYB+2*NGUARD) :: tpdudxc, &
        tpdudyc,tpdvdxc,tpdvdyc


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

  ! write solution data to data.XXXX.XX
  write(filename,'("./IOData/data.",i4.4,".",i6.6,".plt")') count, mype

  i = TecIni('AMR2D'//NULLCHR,'x y u v p denX denY dfun visc curv vort div Q'//NULLCHR,   &
           filename//NULLCHR,'./IOData/'//NULLCHR, &
           Debug,VIsdouble)

  !open(unit=22,file=filename,status='replace')  

  intsx    = (/ (real(i), i=0,NXB) /)
  intsy    = (/ (real(i), i=0,NYB) /)

!  call int2char(mype,index_mype)

  do lb = 1,blockcount


     blockID =  blockList(lb)      


     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)
     dx = del(IAXIS)
     dy = del(JAXIS)
  

     ! Get Coord and Bsize for the block:
     ! Bounding box:
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)

     call Grid_getBlkCenterCoords(blockId,coord)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)


     tpu = 0.
     tpv = 0.
     tpp = 0.
     tpdens = 0.
     tpdensy = 0.

     xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
     xcell = xedge(:) + dx/2.0;
    
     yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
     ycell = yedge(:) + dy/2.0;
    
     facevarxx = facexData(VELC_FACE_VAR,:,:,1)
     facevaryy = faceyData(VELC_FACE_VAR,:,:,1)
 
     facevarr1 = facexData(RH1F_FACE_VAR,:,:,1)
     facevarr2 = faceyData(RH1F_FACE_VAR,:,:,1)
     facevarr3 = facexData(RH2F_FACE_VAR,:,:,1)
     facevarr4 = faceyData(RH2F_FACE_VAR,:,:,1)

     ! U velocity: u(nxb+1,nyb+1)
     ! --------------------------
     tpu = 0.5*(facevarxx(NGUARD+1:nxc,NGUARD:nyc-1)+  &
                facevarxx(NGUARD+1:nxc,NGUARD+1:nyc) )


     ! V velocity: v(nxb+1,nyb+1)
     ! --------------------------                           
     tpv = 0.5*(facevaryy(NGUARD:nxc-1,NGUARD+1:nyc) + &
                facevaryy(NGUARD+1:nxc,NGUARD+1:nyc) )                               

     ! P pressure: p(nxb+1,nyb+1)
     ! -------------------------------
     call centervals2corners(NGUARD,NXB,NYB,nxc,nyc, &
                            solnData(PRES_VAR,:,:,1),tpp)

     ! Distance Function: df(nxb+1,nyb+1)
     ! -------------------------------
     call centervals2corners(NGUARD,NXB,NYB,nxc,nyc, &
                            solnData(DFUN_VAR,:,:,1),tpdfun)

     ! Density: dens(nxb+1,nyb+1)
     ! -------------------------------

     if (dr_nstep .eq. 1) then
        tpdens  = 0.d0
        tpdensy = 0.d0
     else
        tpdens  = 0.5*(1./( facevarr1(NGUARD+1:nxc,NGUARD:nyc-1) + facevarr3(NGUARD+1:nxc,NGUARD:nyc-1) ) +  &
                       1./( facevarr1(NGUARD+1:nxc,NGUARD+1:nyc) + facevarr3(NGUARD+1:nxc,NGUARD+1:nyc) ) )

        tpdensy = 0.5*(1./( facevarr2(NGUARD:nxc-1,NGUARD+1:nyc) + facevarr4(NGUARD:nxc-1,NGUARD+1:nyc) ) +  &
                       1./( facevarr2(NGUARD+1:nxc,NGUARD+1:nyc) + facevarr4(NGUARD+1:nxc,NGUARD+1:nyc) ) )
     end if

     ! Viscosity: visc(nxb+1,nyb+1)
     ! -------------------------------
     call centervals2corners(NGUARD,NXB,NYB,nxc,nyc, &
                            solnData(VISC_VAR,:,:,1),tpvisc)

     ! Viscosity: visc(nxb+1,nyb+1)
     ! -------------------------------
     call centervals2corners(NGUARD,NXB,NYB,nxc,nyc, &
                            solnData(CURV_VAR,:,:,1),tpcurv)

     ! Divergence: 
     ! ----------
     solnData(DUST_VAR,NGUARD:nxc,NGUARD:nyc,1) =      &
             (facevarxx(NGUARD+1:nxc+1,NGUARD:nyc) - &
              facevarxx(NGUARD:nxc,NGUARD:nyc))/dx + &
             (facevaryy(NGUARD:nxc,NGUARD+1:nyc+1) - &
              facevaryy(NGUARD:nxc,NGUARD:nyc))/dy
     call centervals2corners(NGUARD,NXB,NYB,nxc,nyc, &
                             solnData(DUST_VAR,:,:,1),divpp)


            ! Velocity derivatives:
            ! -------- -----------            
!            tpdudxc(ng:nxc,ng:nyc) = (facevarxx(ng+1:nxc+1,ng:nyc) -
!&                                     facevarxx(ng:nxc,ng:nyc))/dx

!            tpdvdyc(ng:nxc,ng:nyc) =  (facevaryy(ng:nxc,ng+1:nyc+1) -
!&                                      facevaryy(ng:nxc,ng:nyc))/dy 

      tpdudycorn(1:NXB+1,1:NYB+1)=(facevarxx(NGUARD+1:nxc,NGUARD+1:nyc)-  &
                                   facevarxx(NGUARD+1:nxc,NGUARD:nyc-1))/dy

      tpdvdxcorn(1:NXB+1,1:NYB+1)=(facevaryy(NGUARD+1:nxc,NGUARD+1:nyc)-  &
                                   facevaryy(NGUARD:nxc-1,NGUARD+1:nyc))/dx 
      
      tpdudxcorn(1:NXB+1,1:NYB+1) = 0.25*( facevarxx(NGUARD+2:nxc+1,NGUARD+1:nyc)+ &
                                           facevarxx(NGUARD+2:nxc+1,NGUARD:nyc-1) - &
                                          facevarxx(NGUARD:nxc-1,NGUARD+1:nyc) - &
                                           facevarxx(NGUARD:nxc-1,NGUARD:nyc-1) )/dx

     tpdvdycorn (1:NXB+1,1:NYB+1) = 0.25*( facevaryy(NGUARD+1:nxc,NGUARD+2:nyc+1) + &
                                           facevaryy(NGUARD:nxc-1,NGUARD+1:nyc+1) - &
                                           facevaryy(NGUARD:nxc-1,NGUARD:nyc-1) - &
                                           facevaryy(NGUARD+1:nxc,NGUARD:nyc-1) )/dy     
   
      ! VORTICITY:
      ! ---------
      ! Corner values of vorticity:
      vortz = tpdvdxcorn - tpdudycorn

     ! Q criterion:
     ! -----------

        qcriterion = -0.5*tpdudxcorn**2 -0.5*tpdvdycorn**2 -tpdudycorn*tpdvdxcorn

     !facevarr1 = facexData(RH1F_FACE_VAR,:,:,1)
     !facevarr2 = facexData(RH2F_FACE_VAR,:,:,1)

     ! vortz = 0.5*( 1. / (facevarr1(NGUARD+1:nxc,NGUARD:nyc-1)  + &
     !                     facevarr2(NGUARD+1:nxc,NGUARD:nyc-1) )+ &
     !               1. / (facevarr1(NGUARD+1:nxc,NGUARD+1:nyc)  + &
!		          facevarr2(NGUARD+1:nxc,NGUARD+1:nyc) ) )

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


      ! Write u:
      arraylb(:,:,1) = sngl(tpu)
      i = TecDat(ijk,arraylb,0)

      ! Write v:
      arraylb(:,:,1) = sngl(tpv)
      i = TecDat(ijk,arraylb,0)

      ! Write p:
      arraylb(:,:,1) = sngl(tpp)
      i = TecDat(ijk,arraylb,0)

      ! Write dens:
      arraylb(:,:,1) = sngl(tpdens)
      i = TecDat(ijk,arraylb,0)

      ! Write dens:
      arraylb(:,:,1) = sngl(tpdensy)
      i = TecDat(ijk,arraylb,0)

      ! Write dfun:
      arraylb(:,:,1) = sngl(tpdfun)
      i = TecDat(ijk,arraylb,0)

      ! Write visc:
      arraylb(:,:,1) = sngl(tpvisc)
      i = TecDat(ijk,arraylb,0)

      ! Write visc:
      arraylb(:,:,1) = sngl(tpcurv)
      i = TecDat(ijk,arraylb,0)


      ! Write omgZ:
      arraylb(:,:,1) = sngl(vortz)
      i = TecDat(ijk,arraylb,0)

      ! Write Div:
      arraylb(:,:,1) = sngl(divpp)
      i = TecDat(ijk,arraylb,0)

      ! Write Q:
      arraylb(:,:,1) = sngl(qcriterion)
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

