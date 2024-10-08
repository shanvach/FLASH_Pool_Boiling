! Subroutine outtotecplot
!
! Subroutine to write out to Tecplot data in binary form.
!
! ---------------------------------------------------------------------------

  subroutine outtotecplot(mype,time,dt,istep,count, &
                          timer,blockList,blockCount,firstfileflag)


  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkPtr, &
    Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_getBlkBoundBox, &
    Grid_getBlkCenterCoords

  use ins_interface, only : ins_velgradtensor

  use Multiphase_data, only: mph_radius

  use Heat_AD_data, only: ht_Nu_l,ht_Nu_t

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
  character(27) :: filename
  character(6) :: index_lb,index_mype

  real xedge(NXB+1),xcell(NXB+1)
  real yedge(NYB+1),ycell(NYB+1)
  real zedge(NZB+1),zcell(NZB+1)
  real intsx(NXB+1), intsy(NYB+1), intsz(NZB+1)

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


  real*4 arraylb(NXB+1,NYB+1,NZB+1)
 
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
  ijk       = (NXB+1)*(NYB+1)*(NZB+1)
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

  if(mype .eq. MASTER_PE) then

   if(firstfileflag .eq. 0) then

        open(unit=44, file="IOData/data_heatFlux",status='replace')
        write(44,55)time,ht_Nu_l,ht_Nu_t,mph_radius
        close(44)

   else

        open(unit=44, file="IOData/data_heatFlux", status='old', position='append')
        write(44,55)time,ht_Nu_l,ht_Nu_t,mph_radius
        close(44)

   end if

  end if

  ! -- data.XXXX.XX --
  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1
  nzc = NZB + NGUARD + 1

! write solution data to data.XXXX.XX

!!-kpd - New Output
!if ( mype .lt. 100) then
!  write(filename,'("./IOData/data.",i4.4,".00",i2.2,".plt")') &
!        count, mype
!elseif ( mype .ge. 100 .AND. mype .lt. 1000) then
!  write(filename,'("./IOData/data.",i4.4,".0",i3.3,".plt")') &
!        count, mype
!else
!  write(filename,'("./IOData/data.",i4.4,".",i4.4,".plt")') &
!        count, mype
!end if

  write(filename,'("./IOData/data.",i4.4,".",i4.4,".plt")') &
        count, mype

!  write(filename,'("./IOData/data.",i4.4,".",i3.3,".plt")') &
!        count, mype


  i = TecIni('AMR3D'//NULLCHR,                            &
             'x y z source analytical numerical error coeffx coeffy coeffz'//NULLCHR,  &
             filename//NULLCHR,                           &
             './IOData/'//NULLCHR,                  &
             Debug,VIsdouble)



  ! open(unit=22,file=filename,status='replace')  

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

     xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
     xcell = xedge(:) + dx/2.0;
    
     yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
     ycell = yedge(:) + dy/2.0;
    
     zedge = coord(KAXIS) - bsize(KAXIS)/2.0 + dz*intsz;
     zcell = zedge(:) + dz/2.0; 


     facevarxx = facexData(RH1F_FACE_VAR,:,:,:) + facexData(RH2F_FACE_VAR,:,:,:)
     facevaryy = faceyData(RH1F_FACE_VAR,:,:,:) + faceyData(RH2F_FACE_VAR,:,:,:)
     facevarzz = facezData(RH1F_FACE_VAR,:,:,:) + facezData(RH2F_FACE_VAR,:,:,:)         

     ! U velocity: u(NXB+1,NYB+1,NZB+1)
     ! --------------------------------
     vortx = 0.25*(facevarxx(NGUARD+1:nxc,NGUARD:nyc-1,NGUARD:nzc-1)+ &
                 facevarxx(NGUARD+1:nxc,NGUARD:nyc-1,NGUARD+1:nzc)+ &
                 facevarxx(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc)+ &
                 facevarxx(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD:nzc-1))
    
     ! V velocity: v(NXB+1,NYB+1,NZB+1)
     ! --------------------------------                           
     vorty = 0.25*(facevaryy(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD:nzc-1) + &
                 facevaryy(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD:nzc-1) + &
                 facevaryy(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc) + &
                 facevaryy(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD+1:nzc));
    
     vortz = 0.25*(facevarzz(NGUARD:nxc-1,NGUARD:nyc-1,NGUARD+1:nzc)   + &
                 facevarzz(NGUARD+1:nxc,NGUARD:nyc-1,NGUARD+1:nzc)   + &
                 facevarzz(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc)   + &
                 facevarzz(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD+1:nzc));
      

     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc, &
                             solnData(PRES_VAR,:,:,:),tpp)


     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             solnData(RTES_VAR,:,:,:),TVtpp)


     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             solnData(DELP_VAR,:,:,:),tempp)

     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             solnData(PTES_VAR,:,:,:),mdotp)

    ! Write Block Results into data file:
     call int2char(lb,index_lb)

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

     ! Write TV:
     arraylb = sngl(TVtpp)
     i = TecDat(ijk,arraylb,0)

     ! Write p:
     arraylb = sngl(tpp)
     i = TecDat(ijk,arraylb,0)

     arraylb = sngl(tempp)
     i = TecDat(ijk,arraylb,0)

     arraylb = sngl(mdotp)
     i = TecDat(ijk,arraylb,0)

     arraylb = sngl(vortx)
     i = TecDat(ijk,arraylb,0)

     arraylb = sngl(vorty)
     i = TecDat(ijk,arraylb,0)

     arraylb = sngl(vortz)
     i = TecDat(ijk,arraylb,0)

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
55    format(g23.15,g23.15,g23.15,g23.15)

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


