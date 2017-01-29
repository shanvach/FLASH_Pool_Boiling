! Subroutine outtotecplot
!
! Subroutine to write out to Tecplot data in binary form.
!
! ---------------------------------------------------------------------------
#include "constants.h"
#include "Flash.h"

  subroutine outtotecplot(mype,time,dt,istep,count, &
                          timer,blockList,blockCount,firstfileflag)


  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkPtr, &
                 Grid_releaseBlkPtr, Grid_getBlkIndexLimits, &
                 Grid_getBlkBoundBox,Grid_getBlkCenterCoords

  use ins_interface, only : ins_velgradtensor

#ifdef FLASH_GRID_PARAMESH
  use physicaldata, only : interp_mask_unk, interp_mask_unk_res
#endif
  implicit none
#include "Flash_mpi.h"
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
            Sxy,Syz,Sxz,Oxy,Oyz,Oxz,Qcr,divpp,TVtpp


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

  logical :: file_exists

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
  INQUIRE(FILE=filename, EXIST=file_exists)

  ! create/clear filetime.XX if time = 0
  if( (firstfileflag .eq. 0) .or. (.not. file_exists) ) then
     open(unit=33, file=filename, status='unknown')

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

! write solution data to data.XXXX.XX

  write(filename,'("./IOData/data.",i4.4,".",i4.4,".plt")') &
        count, mype


  i = TecIni('AMR3D'//NULLCHR,                            &
             'x y z u v w p wx wy wz Q div TV'//NULLCHR,  &
             filename//NULLCHR,                           &
             './IOData/'//NULLCHR,                  &
             Debug,VIsdouble)




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

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     tpu = 0.
     tpv = 0.
     tpw = 0.
     tpp = 0.
     omg = 0.

     xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
     xcell = xedge(:) + dx/2.0;
    
     yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
     ycell = yedge(:) + dy/2.0;
    
     zedge = coord(KAXIS) - bsize(KAXIS)/2.0 + dz*intsz;
     zcell = zedge(:) + dz/2.0; 


     facevarxx = facexData(VELC_FACE_VAR,:,:,:)
     facevaryy = faceyData(VELC_FACE_VAR,:,:,:)
     facevarzz = facezData(VELC_FACE_VAR,:,:,:)          
 
     ! U velocity: u(NXB+1,NYB+1,NZB+1)
     ! --------------------------------
     tpu = 0.25*(facevarxx(NGUARD+1:nxc,NGUARD:nyc-1,NGUARD:nzc-1)+ &
                 facevarxx(NGUARD+1:nxc,NGUARD:nyc-1,NGUARD+1:nzc)+ &
                 facevarxx(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc)+ &
                 facevarxx(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD:nzc-1))
    
!!$     ! y = 1, z = 1:    
!!$     tpu(:,1,1) = .5*(facevarxx(NGUARD+1:nxc,NGUARD,NGUARD+1) + &
!!$                      facevarxx(NGUARD+1:nxc,NGUARD+1,NGUARD));
!!$                              
!!$     ! y = 1, z = end:    
!!$     tpu(:,1,NZB+1) = .5*(facevarxx(NGUARD+1:nxc,NGUARD+1,nzc)+ &
!!$                          facevarxx(NGUARD+1:nxc,NGUARD,nzc-1)) 
!!$                              
!!$     ! y = end, z = 1:    
!!$     tpu(:,NYB+1,1) = .5*(facevarxx(NGUARD+1:nxc,nyc,NGUARD+1)+ &
!!$                          facevarxx(NGUARD+1:nxc,nyc-1,NGUARD))        
!!$
!!$                              
!!$     ! y = end, z = end:    
!!$     tpu(:,NYB+1,NZB+1)=.5*(facevarxx(NGUARD+1:nxc,nyc,nzc-1)+ &
!!$                            facevarxx(NGUARD+1:nxc,nyc-1,nzc)) 


     ! V velocity: v(NXB+1,NYB+1,NZB+1)
     ! --------------------------------                           
     tpv = 0.25*(facevaryy(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD:nzc-1) + &
                 facevaryy(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD:nzc-1) + &
                 facevaryy(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc) + &
                 facevaryy(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD+1:nzc));
    
!!$     !x = 1, z = 1:    
!!$     tpv(1,:,1) = .5*(facevaryy(NGUARD,NGUARD+1:nyc,NGUARD+1) + &
!!$                      facevaryy(NGUARD+1,NGUARD+1:nyc,NGUARD));
!!$                              
!!$     ! x = 1, z = end:    
!!$     tpv(1,:,NZB+1) = .5*(facevaryy(NGUARD+1,NGUARD+1:nyc,nzc) + &
!!$                          facevaryy(NGUARD,NGUARD+1:nyc,nzc-1)); 
!!$                              
!!$     ! x = end, z = 1:    
!!$     tpv(NXB+1,:,1) = .5*(facevaryy(nxc,NGUARD+1:nyc,NGUARD+1) + &
!!$                          facevaryy(nxc-1,NGUARD+1:nyc,NGUARD));        
!!$                              
!!$     ! x = end, z = end:    
!!$     tpv(NXB+1,:,NZB+1) = .5*(facevaryy(nxc,NGUARD+1:nyc,nzc-1) + &
!!$                              facevaryy(nxc-1,NGUARD+1:nyc,nzc));         

     ! W velocity: w(NXB+1,NYB+1,NZB+1)
     ! --------------------------------
     tpw = 0.25*(facevarzz(NGUARD:nxc-1,NGUARD:nyc-1,NGUARD+1:nzc)   + &
                 facevarzz(NGUARD+1:nxc,NGUARD:nyc-1,NGUARD+1:nzc)   + &
                 facevarzz(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc)   + &
                 facevarzz(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD+1:nzc));
      
        
!!$     ! x = 1, y = 1:    
!!$     tpw(1,1,:) = .5*(facevarzz(NGUARD,NGUARD+1,NGUARD+1:nzc) + &
!!$                      facevarzz(NGUARD+1,NGUARD,NGUARD+1:nzc));
!!$                              
!!$     ! x = 1, y = end:    
!!$     tpw(1,NYB+1,:) = .5*(facevarzz(NGUARD+1,nyc,NGUARD+1:nzc) + &
!!$                          facevarzz(NGUARD,nyc-1,NGUARD+1:nzc)); 
!!$                              
!!$     ! x = end, y = 1:    
!!$     tpw(NXB+1,1,:) = .5*(facevarzz(nxc,NGUARD+1,NGUARD+1:nzc) + &
!!$                          facevarzz(nxc-1,NGUARD,NGUARD+1:nzc));        
!!$
!!$                              
!!$     ! x = end, y = end:    
!!$     tpw(NXB+1,NYB+1,:) = .5*(facevarzz(nxc,nyc-1,NGUARD+1:nzc) + &
!!$                              facevarzz(nxc-1,nyc,NGUARD+1:nzc));       
                              

     ! P pressure: p(NXB+1,NYB+1,NZB+1)
     ! -------------------------------
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc, &
                             solnData(PRES_VAR,:,:,:),tpp)


     ! Divergence: ! Ojo Intermediate velocities: unk(3, for div(u) unk(5 !!!!!!!
     ! ----------
!     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
!                             unk(3,:,:,:,lb),divppstar)

     ! TV : TV(NXB+1,NYB+1,NZB+1)
     ! -------------------------------
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             solnData(TVIS_VAR,:,:,:),TVtpp)


     ! Divergence: 
     ! ----------
     solnData(DUST_VAR,:,:,:) = 0.
     solnData(DUST_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1) =      &
             (facevarxx(NGUARD+2:nxc  ,NGUARD+1:nyc-1,NGUARD+1:nzc-1)     - &
              facevarxx(NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1))/dx + &
             (facevaryy(NGUARD+1:nxc-1,NGUARD+2:nyc  ,NGUARD+1:nzc-1)     - &
              facevaryy(NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1))/dy + &
             (facevarzz(NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+2:nzc  )     - &
              facevarzz(NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1))/dz 
    

     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             solnData(DUST_VAR,:,:,:),divpp)


     ! Velocity derivatives:
     ! -------- -----------            
     call ins_velgradtensor(NGUARD,facexData,faceyData,facezData, &
                               dx,dy,dz,tpdudxc,tpdudyc,tpdudzc,&
                tpdvdxc,tpdvdyc,tpdvdzc,tpdwdxc,tpdwdyc,tpdwdzc)

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

!!$     do k = 1,NZB+1
!!$        do j = 1,NYB+1
!!$           do i = 1,NXB+1
!!$              omg(i,j,k) = sqrt(vortx(i,j,k)**2 +  &
!!$                                vorty(i,j,k)**2 +  &
!!$                                vortz(i,j,k)**2)
!!$           enddo
!!$        enddo
!!$     enddo

     ! Q criterion:
     ! First Sij:
     !Sxx = tpdudxcorn
     Sxy = 0.5*(tpdudycorn + tpdvdxcorn)
     Sxz = 0.5*(tpdudzcorn + tpdwdxcorn)
     !Syy = tpdvdycorn
     Syz = 0.5*(tpdvdzcorn + tpdwdycorn)
     !Szz = tpdwdzcorn

     ! Then Oij:
     Oxy = 0.5*(tpdudycorn - tpdvdxcorn)
     Oxz = 0.5*(tpdudzcorn - tpdwdxcorn)
     Oyz = 0.5*(tpdvdzcorn - tpdwdycorn)

     ! Calculate Q:
     Qcr = 0.5*(2.*(Oxy*Oxy+Oyz*Oyz+Oxz*Oxz) -                     &
               (tpdudxcorn*tpdudxcorn + tpdvdycorn*tpdvdycorn +    &
                tpdwdzcorn*tpdwdzcorn +                            &
                2.*(Sxy*Sxy+Syz*Syz+Sxz*Sxz)))                      


     ! Write Block Results into data file:
     write(index_lb,"(I6.6)") blockID
     i = TecZne(                                                       &
                'ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR,  &
                 NXB+1,NYB+1,NZB+1,                                    &
                 'BLOCK'//NULLCHR,                                     &
                 CHAR(0))

            
     ! Write x:
     do k=1,NZB+1
        do j=1,NYB+1
           do i=1,NXB+1
              arraylb(i,j,k) = real(xedge(i),KIND=4) !sngl(xedge(i))
           enddo
        enddo
     enddo
     i = TecDat(ijk,arraylb,0)


     ! Write y:
     do k=1,NZB+1
        do j=1,NYB+1
           do i=1,NXB+1
              arraylb(i,j,k) = real(yedge(j),KIND=4) !sngl(yedge(j))
           enddo
        enddo
     enddo
     i = TecDat(ijk,arraylb,0)


     ! Write z:
     do k=1,NZB+1
        do j=1,NYB+1
           do i=1,NXB+1
              arraylb(i,j,k) = real(zedge(k),KIND=4) !sngl(zedge(k))
           enddo
        enddo
     enddo
     i = TecDat(ijk,arraylb,0)


     ! Write u:
     arraylb = real(tpu,KIND=4) !sngl(tpu)
     i = TecDat(ijk,arraylb,0)

     ! Write v:
     arraylb = real(tpv,KIND=4) !sngl(tpv)
     i = TecDat(ijk,arraylb,0)

     ! Write w:
     arraylb = real(tpw,KIND=4) !sngl(tpw)
     i = TecDat(ijk,arraylb,0)

     ! Write p:
     arraylb = real(tpp,KIND=4) !sngl(tpp)
     i = TecDat(ijk,arraylb,0)

     ! Write omgX:
     arraylb = real(vortx,KIND=4) !sngl(vortx)
     i = TecDat(ijk,arraylb,0)


     ! Write omgY:
     arraylb = real(vorty,KIND=4) !sngl(vorty)
     i = TecDat(ijk,arraylb,0)


     ! Write omgZ:
     arraylb = real(vortz,KIND=4) !sngl(vortz)
     i = TecDat(ijk,arraylb,0)

     ! Write omg:
!!$     do k=1,NZB+1
!!$        do j=1,NYB+1
!!$           do i=1,NXB+1
!!$              write(22,'g14.8') omg(i,j,k)
!!$           enddo
!!$        enddo
!!$     enddo

     ! Write Q:
     arraylb = real(Qcr,KIND=4) !sngl(Qcr)
     i = TecDat(ijk,arraylb,0)

     ! Write Div:
     arraylb = real(divpp,KIND=4) !sngl(divpp)
     i = TecDat(ijk,arraylb,0)

     ! Write TV:
     arraylb = real(TVtpp,KIND=4) !sngl(TVtpp)
     i = TecDat(ijk,arraylb,0)

     ! Write Div Ustar:
!!$     do k=1,NZB+1
!!$        do j=1,NYB+1
!!$           do i=1,NXB+1
!!$              write(22,'g14.8') divppstar(i,j,k)
!!$           enddo
!!$        enddo
!!$     enddo


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

