! Subroutine outtotecplot
!
! Subroutine to write out to Tecplot data in binary form.
!
! ---------------------------------------------------------------------------

#include "constants.h"
#include "Flash.h"


  subroutine outtoflowvc(mype,time,dt,istep,count,&
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
           vortz,divpp,tpdens,tpdensy,tpdfun,tpvisc,tpcurv

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


  ! -- data.XXXX.XX --
  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1

  ! write solution data to data.XXXX.XX
  write(filename,'("./IOData/data.",i4.4,".",i6.6,".bin")') count, mype

!  i = TecIni('AMR2D'//NULLCHR,'x y u v p denX denY dfun visc curv vort div'//NULLCHR,   &
!           filename//NULLCHR,'./IOData/'//NULLCHR, &
!           Debug,VIsdouble)

  open(unit=22,file=filename,status='replace')  
  close(22)

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
         
      ! VORTICITY:
      ! ---------
      ! Corner values of vorticity:
      vortz = tpdvdxcorn - tpdudycorn

     !facevarr1 = facexData(RH1F_FACE_VAR,:,:,1)
     !facevarr2 = facexData(RH2F_FACE_VAR,:,:,1)

     ! vortz = 0.5*( 1. / (facevarr1(NGUARD+1:nxc,NGUARD:nyc-1)  + &
     !                     facevarr2(NGUARD+1:nxc,NGUARD:nyc-1) )+ &
     !               1. / (facevarr1(NGUARD+1:nxc,NGUARD+1:nyc)  + &
!		          facevarr2(NGUARD+1:nxc,NGUARD+1:nyc) ) )

      ! Write Block Results into data file:
      call int2char(lb,index_lb)

      !i = TecZne('ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR, &
      !    NXB+1,NYB+1,1,'BLOCK'//NULLCHR,CHAR(0))

            
      ! Write x:
      do j=1,NYB+1
         do i=1,NXB+1
            !arraylb(i,j,1) = sngl(xedge(i))
         enddo
      enddo
      !i = TecDat(ijk,arraylb,0)


      ! Write y:
      do j=1,NYB+1
         do i=1,NXB+1
            !arraylb(i,j,1) = sngl(yedge(j))
         enddo
      enddo
      !i = TecDat(ijk,arraylb,0)


      ! Write u:
      !arraylb(:,:,1) = sngl(tpu)
      !i = TecDat(ijk,arraylb,0)

      ! Write v:
      !arraylb(:,:,1) = sngl(tpv)
      !i = TecDat(ijk,arraylb,0)

      ! Write p:
      !arraylb(:,:,1) = sngl(tpp)
      !i = TecDat(ijk,arraylb,0)

      ! Write dens:
      !arraylb(:,:,1) = sngl(tpdens)
      !i = TecDat(ijk,arraylb,0)

      ! Write dens:
      !arraylb(:,:,1) = sngl(tpdensy)
      !i = TecDat(ijk,arraylb,0)

      ! Write dfun:
      !arraylb(:,:,1) = sngl(tpdfun)
      !i = TecDat(ijk,arraylb,0)

      ! Write visc:
      !arraylb(:,:,1) = sngl(tpvisc)
      !i = TecDat(ijk,arraylb,0)

      ! Write visc:
      !arraylb(:,:,1) = sngl(tpcurv)
      !i = TecDat(ijk,arraylb,0)


      ! Write omgZ:
      !arraylb(:,:,1) = sngl(vortz)
      !i = TecDat(ijk,arraylb,0)

      ! Write Div:
      !arraylb(:,:,1) = sngl(divpp)
      !i = TecDat(ijk,arraylb,0)

   enddo

   !i = TecEnd()


  End subroutine outtoflowvc

