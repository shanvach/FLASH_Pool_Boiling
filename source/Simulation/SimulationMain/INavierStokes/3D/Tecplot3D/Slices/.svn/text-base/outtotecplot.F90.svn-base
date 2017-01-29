! Subroutine outtotecplot
!
! Subroutine to write out to Tecplot data in binary form.
!
! ---------------------------------------------------------------------------
#include "constants.h"
#include "Flash.h"

#define ALLOC_TECPLOT_VARS
!!#define ONE_SLC_FILE

  subroutine outtotecplot(mype,time,dt,istep,count, &
                          timer,blockList,blockCount,firstfileflag)


  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkPtr, &
                 Grid_releaseBlkPtr, Grid_getBlkIndexLimits, &
                 Grid_getBlkBoundBox,Grid_getBlkCenterCoords

  use ins_interface, only : ins_velgradtensor

  use Grid_data, only : gr_meshComm,gr_meshnumProcs

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

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
  character(25) :: filename
  character(32) :: filename2
  character(6) :: index_lb,index_mype

  real xedge(NXB+1),xcell(NXB+1)
  real yedge(NYB+1),ycell(NYB+1)
  real zedge(NZB+1),zcell(NZB+1)
  real intsx(NXB+1), intsy(NYB+1), intsz(NZB+1)

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData


#ifdef ALLOC_TECPLOT_VARS

  real, allocatable, dimension(:,:,:) :: facevarxx,facevaryy,facevarzz,tpu,tpv,tpw,tpp, &
            tpdudxcorn, tpdudycorn, tpdudzcorn,&
            tpdvdxcorn, tpdvdycorn, tpdvdzcorn,&
            tpdwdxcorn, tpdwdycorn, tpdwdzcorn,&
            vortx,vorty,vortz,omg,             &
            Sxy,Syz,Sxz,Oxy,Oyz,Oxz,Qcr,divpp, &
            TVtpp, tpdudxc, tpdudyc,tpdudzc,   &
            tpdvdxc,tpdvdyc,tpdvdzc,tpdwdxc,   &
            tpdwdyc,tpdwdzc

  real*4, allocatable, dimension(:,:,:) :: arraylb 

#else

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

#endif

  integer blockID

  real del(3),dx,dy,dz
  real, dimension(MDIM)  :: coord,bsize
  real ::  boundBox(2,MDIM)

  real, save :: tp_xmin,tp_xmax,tp_ymin,tp_ymax,tp_zmin,tp_zmax

  integer, parameter :: XPOS_TP = 1
  integer, parameter :: YPOS_TP = 2
  integer, parameter :: ZPOS_TP = 3
  integer, parameter :: UVEL_TP = 4
  integer, parameter :: VVEL_TP = 5
  integer, parameter :: WVEL_TP = 6
  integer, parameter :: PRES_TP = 7
  integer, parameter :: XVOR_TP = 8  
  integer, parameter :: YVOR_TP = 9  
  integer, parameter :: ZVOR_TP =10
  integer, parameter :: QCRT_TP =12
  integer, parameter :: NUM_EXTP=12      


  integer, save :: nslicex
  real, save, allocatable, dimension(:) :: Xslice 
  integer, save :: nslicey
  real, save, allocatable, dimension(:) :: Yslice 
  integer, save :: nslicez
  real, save, allocatable, dimension(:) :: Zslice 

  logical, save :: readfilex,readfiley,readfilez


  integer :: icellx,islicex,slcx_zone,slcx_zone_total
  integer :: icelly,islicey,slcy_zone,slcy_zone_total
  integer :: icellz,islicez,slcz_zone,slcz_zone_total  
  real :: Xslci,Yslci,Zslci,intFactor_low,intFactor_high

  real*4, allocatable, dimension(:,:,:,:) :: slicedata  

  integer :: iproc,izone,buffer_size,ierr 

  integer,dimension(MPI_STATUS_SIZE) :: status

  integer*4 TecIni,TecDat,TecZne,TecNod,TecFil,TecEnd
  integer*4 VIsdouble
  integer*4 Debug,ijkx,ijky,ijkz,Npts,NElm
  character*1 NULLCHR

  logical, save :: firstcall=.true.

!-----------------------------------------------------------------------
!                                                         TecPlot set-up
!-----------------------------------------------------------------------
  Debug     = 0
  VIsdouble = 0
  NULLCHR   = CHAR(0)
  ijkx      = (1)*(NYB+1)*(NZB+1)
  ijky      = (NXB+1)*(1)*(NZB+1)
  ijkz      = (NXB+1)*(NYB+1)*(1)
!-----------------------------------------------------------------------

#ifdef ALLOC_TECPLOT_VARS

  allocate(facevarxx(NXB+2*NGUARD+1,NYB+2*NGUARD,NZB+2*NGUARD))
  allocate(facevaryy(NXB+2*NGUARD,NYB+2*NGUARD+1,NZB+2*NGUARD))
  allocate(facevarzz(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD+1))

  allocate(tpu(NXB+1,NYB+1,NZB+1),tpv(NXB+1,NYB+1,NZB+1),tpw(NXB+1,NYB+1,NZB+1))
  allocate(tpp(NXB+1,NYB+1,NZB+1))
  allocate(tpdudxcorn(NXB+1,NYB+1,NZB+1),tpdudycorn(NXB+1,NYB+1,NZB+1),tpdudzcorn(NXB+1,NYB+1,NZB+1))
  allocate(tpdvdxcorn(NXB+1,NYB+1,NZB+1),tpdvdycorn(NXB+1,NYB+1,NZB+1),tpdvdzcorn(NXB+1,NYB+1,NZB+1))
  allocate(tpdwdxcorn(NXB+1,NYB+1,NZB+1),tpdwdycorn(NXB+1,NYB+1,NZB+1),tpdwdzcorn(NXB+1,NYB+1,NZB+1))
  allocate(vortx(NXB+1,NYB+1,NZB+1),vorty(NXB+1,NYB+1,NZB+1),vortz(NXB+1,NYB+1,NZB+1))
  allocate(Sxy(NXB+1,NYB+1,NZB+1),Syz(NXB+1,NYB+1,NZB+1),Sxz(NXB+1,NYB+1,NZB+1))
  allocate(Oxy(NXB+1,NYB+1,NZB+1),Oyz(NXB+1,NYB+1,NZB+1),Oxz(NXB+1,NYB+1,NZB+1),Qcr(NXB+1,NYB+1,NZB+1))

  allocate(omg(NXB+1,NYB+1,NZB+1),divpp(NXB+1,NYB+1,NZB+1),TVtpp(NXB+1,NYB+1,NZB+1))

  allocate(arraylb(NXB+1,NYB+1,NZB+1))

  allocate(tpdudxc(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD), & 
           tpdudyc(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD), &
           tpdudzc(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD))

  allocate(tpdvdxc(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD), &
           tpdvdyc(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD), &
           tpdvdzc(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD))

  allocate(tpdwdxc(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD), &
           tpdwdyc(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD), &
           tpdwdzc(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD))
#endif

  if (firstcall) then

     call RuntimeParameters_get("xmin", tp_xmin)
     call RuntimeParameters_get("xmax", tp_xmax)     
     call RuntimeParameters_get("ymin", tp_ymin)
     call RuntimeParameters_get("ymax", tp_ymax)
     call RuntimeParameters_get("zmin", tp_zmin)
     call RuntimeParameters_get("zmax", tp_zmax)

     ! Here we read slices:

     ! X slices:
     call RuntimeParameters_get("tecplot_nslicesx", nslicex)     

     if (nslicex .gt. 0) then
        allocate(Xslice(nslicex))
        call RuntimeParameters_get("tecplot_slicesx_from_file", readfilex)      
        if (readfilex) then
           open(unit=33, file='./xlist.slc', status='old')
           do islicex = 1,nslicex
              read(33,*) Xslice(islicex)
           end do
           close(33)
        else
           ! write evenly spaced across the domain:
           dx = (tp_xmax-tp_xmin)/real(nslicex)
           do islicex =1,nslicex
              Xslice(islicex) = tp_xmin + 0.5*dx + real(islicex-1)*dx
           enddo
        endif
     endif


     ! Y slices:
     call RuntimeParameters_get("tecplot_nslicesy", nslicey)     

     if (nslicey .gt. 0) then
        allocate(Yslice(nslicey))
        call RuntimeParameters_get("tecplot_slicesy_from_file", readfiley)      
        if (readfiley) then
           open(unit=33, file='./ylist.slc', status='old')
           do islicey = 1,nslicey
              read(33,*) Yslice(islicey)
           end do
           close(33)
        else
           ! write evenly spaced across the domain:
           dy = (tp_ymax-tp_ymin)/real(nslicey)
           do islicey =1,nslicey
              Yslice(islicey) = tp_ymin + 0.5*dy + real(islicey-1)*dy
           enddo
        endif
     endif


     ! Z slices:
     call RuntimeParameters_get("tecplot_nslicesz", nslicez)

     if (nslicez .gt. 0) then
        allocate(Zslice(nslicez))
        call RuntimeParameters_get("tecplot_slicesz_from_file", readfilez)      
        if (readfilez) then
           open(unit=33, file='./zlist.slc', status='old')
           do islicez = 1,nslicez
              read(33,*) Zslice(islicez)
           end do
           close(33)
        else
           ! write evenly spaced across the domain:
           dz = (tp_zmax-tp_zmin)/real(nslicez)
           do islicez =1,nslicez
              Zslice(islicez) = tp_zmin + 0.5*dz + real(islicez-1)*dz
           enddo
        endif
     endif

     firstcall = .false.

  endif

#ifdef FLASH_GRID_UG
  if (mype .eq. MASTER_PE) then
#endif
! -- filetime.XX --
  write(filename, '("IOData/data_time.", i5.5)') mype
  ! create/clear filetime.XX if time = 0
  if(firstfileflag .eq. 0) then
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
#ifdef FLASH_GRID_UG
  endif
#endif

  ! -- data.XXXX.XX --
  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1
  nzc = NZB + NGUARD + 1

  ! Now all processors compute and interpolate to the slice:
  intsx    = (/ (real(i), i=0,NXB) /)
  intsy    = (/ (real(i), i=0,NYB) /)
  intsz    = (/ (real(i), i=0,NZB) /)

#ifndef FLASH_GRID_PARAMESH

  lb      = CONSTANT_ONE
  blockID = blockList(lb)

  ! Get Coord and Bsize for the block:
  ! Bounding box:
  call Grid_getBlkBoundBox(blockId,boundBox)
  
  ! Get blocks dx, dy ,dz:
  call Grid_getDeltas(blockID,del)
  dx = del(IAXIS)
  dy = del(JAXIS)
  dz = del(KAXIS)
  
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

  ! V velocity: v(NXB+1,NYB+1,NZB+1)
  ! --------------------------------                           
  tpv = 0.25*(facevaryy(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD:nzc-1) + &
              facevaryy(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD:nzc-1) + &
              facevaryy(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc) + &
              facevaryy(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD+1:nzc));      

  ! W velocity: w(NXB+1,NYB+1,NZB+1)
  ! --------------------------------
  tpw = 0.25*(facevarzz(NGUARD:nxc-1,NGUARD:nyc-1,NGUARD+1:nzc)   + &
              facevarzz(NGUARD+1:nxc,NGUARD:nyc-1,NGUARD+1:nzc)   + &
              facevarzz(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc)   + &
              facevarzz(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD+1:nzc)); 

  ! P pressure: p(NXB+1,NYB+1,NZB+1)
  ! -------------------------------
  call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc, &
                          solnData(PRES_VAR,:,:,:),tpp)

  ! Divergence: ! Ojo Intermediate velocities: unk(3, for div(u) unk(5 !!!!!!!
  ! ----------
!  call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
!                          unk(3,:,:,:,lb),divppstar)

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
  call ins_velgradtensor(NGUARD,facexData,faceyData,facezData,           &
                         dx,dy,dz,tpdudxc,tpdudyc,tpdudzc,               &
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
  Qcr = 0.5*(2.*(Oxy*Oxy+Oyz*Oyz+Oxz*Oxz) -                &
       (tpdudxcorn*tpdudxcorn + tpdvdycorn*tpdvdycorn +    &
       tpdwdzcorn*tpdwdzcorn + 2.*(Sxy*Sxy+Syz*Syz+Sxz*Sxz)))  

  ! Release Pointers:
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
  call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

#endif

  write(index_mype,"(I6.6)") mype


  ! Write solution to x slices:
  ! write solution data to slcx.XXXX.XX.plt
  allocate(slicedata(NYB+1,NZB+1,NUM_EXTP,MAXBLOCKS)); slicedata = 0.
  do islicex =1,nslicex

     Xslci = Xslice(islicex) 

#ifdef ONE_SLC_FILE
     if (mype .eq. MASTER_PE) then
        write(filename,'("./IOData/slcx.",i4.4,".",i2.2,".plt")') &
              count, islicex

        i = TecIni('AMR3D'//NULLCHR,                            &
                   'x y z u v w p wx wy wz Q'//NULLCHR,         & !div TV
                   filename//NULLCHR,                           &
                   './IOData/'//NULLCHR,                        &
                   Debug,VIsdouble)
     endif
#endif


     ! Loop over blocks:
     slcx_zone = 0
     do lb = 1,blockcount
  
        blockID =  blockList(lb)      

#ifdef FLASH_GRID_PARAMESH
        ! Get Coord and Bsize for the block:
        ! Bounding box:
        call Grid_getBlkBoundBox(blockId,boundBox)
#endif

        ! Test if slice falls on block domain
        if (Xslci .gt. boundBox(2,IAXIS)) cycle
        if (Xslci .le. boundBox(1,IAXIS)) cycle 

#ifdef FLASH_GRID_PARAMESH
        ! Get blocks dx, dy ,dz:
        call Grid_getDeltas(blockID,del)
        dx = del(IAXIS)
        dy = del(JAXIS)
        dz = del(KAXIS)
        
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

        ! V velocity: v(NXB+1,NYB+1,NZB+1)
        ! --------------------------------                           
        tpv = 0.25*(facevaryy(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD:nzc-1) + &
             facevaryy(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD:nzc-1) + &
             facevaryy(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc) + &
             facevaryy(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD+1:nzc));      

        ! W velocity: w(NXB+1,NYB+1,NZB+1)
        ! --------------------------------
        tpw = 0.25*(facevarzz(NGUARD:nxc-1,NGUARD:nyc-1,NGUARD+1:nzc)   + &
             facevarzz(NGUARD+1:nxc,NGUARD:nyc-1,NGUARD+1:nzc)   + &
             facevarzz(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc)   + &
             facevarzz(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD+1:nzc)); 

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
        call ins_velgradtensor(NGUARD,facexData,faceyData,facezData,           &
                               dx,dy,dz,tpdudxc,tpdudyc,tpdudzc,               &
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
        Qcr = 0.5*(2.*(Oxy*Oxy+Oyz*Oyz+Oxz*Oxz) -                &
             (tpdudxcorn*tpdudxcorn + tpdvdycorn*tpdvdycorn +    &
              tpdwdzcorn*tpdwdzcorn + 2.*(Sxy*Sxy+Syz*Syz+Sxz*Sxz)))                      

        ! Release Pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

        ! Interpolate values to slice location:
        icellx = floor((Xslci - xedge(1))/dx) + 1
        if (icellx .eq. (NXB+1)) icellx = icellx - 1

        ! Linear interpolation factors
        intFactor_low  = (xedge(icellx+1)-Xslci)/dx
        intFactor_high = (Xslci - xedge(icellx))/dx

        ! Interpolate to slice array:
        slcx_zone = slcx_zone + 1
        ! X position:
        do j=1,NZB+1
           do i=1,NYB+1
              slicedata(i,j,XPOS_TP,slcx_zone)= &
                   real(intFactor_low*xedge(icellx)+intFactor_high*xedge(icellx+1),KIND=4)
           enddo
        enddo
        
        ! Y Position:
        do j=1,NZB+1
           do i=1,NYB+1
              slicedata(i,j,YPOS_TP,slcx_zone)=real(yedge(i),KIND=4)
           enddo
        enddo
        
        ! Z Position:
        do j=1,NZB+1
           do i=1,NYB+1
              slicedata(i,j,ZPOS_TP,slcx_zone)=real(zedge(j),KIND=4)
           enddo
        enddo
        
        ! U velocity:
        do j=1,NZB+1
           do i=1,NYB+1
              slicedata(i,j,UVEL_TP,slcx_zone)= &
                   real(intFactor_low*tpu(icellx,i,j)+intFactor_high*tpu(icellx+1,i,j),KIND=4)
           enddo
        enddo
        
        ! V velocity:
        do j=1,NZB+1
           do i=1,NYB+1
              slicedata(i,j,VVEL_TP,slcx_zone)= &
                   real(intFactor_low*tpv(icellx,i,j)+intFactor_high*tpv(icellx+1,i,j),KIND=4)
           enddo
        enddo
        
        ! W velocity:
        do j=1,NZB+1
           do i=1,NYB+1
              slicedata(i,j,WVEL_TP,slcx_zone)= &
                   real(intFactor_low*tpw(icellx,i,j)+intFactor_high*tpw(icellx+1,i,j),KIND=4)
           enddo
        enddo
        
        ! Pressure:
        do j=1,NZB+1
           do i=1,NYB+1
              slicedata(i,j,PRES_TP,slcx_zone)= &
                   real(intFactor_low*tpp(icellx,i,j)+intFactor_high*tpp(icellx+1,i,j),KIND=4)
           enddo
        enddo
        
        ! Vortx:
        do j=1,NZB+1
           do i=1,NYB+1
              slicedata(i,j,XVOR_TP,slcx_zone)= &
                   real(intFactor_low*vortx(icellx,i,j)+intFactor_high*vortx(icellx+1,i,j),KIND=4)
           enddo
        enddo
        
        ! Vorty:
        do j=1,NZB+1
           do i=1,NYB+1
              slicedata(i,j,YVOR_TP,slcx_zone)= &
                   real(intFactor_low*vorty(icellx,i,j)+intFactor_high*vorty(icellx+1,i,j),KIND=4)
           enddo
        enddo
        
        ! Vortz:
        do j=1,NZB+1
           do i=1,NYB+1
              slicedata(i,j,ZVOR_TP,slcx_zone)= &
                   real(intFactor_low*vortz(icellx,i,j)+intFactor_high*vortz(icellx+1,i,j),KIND=4)
           enddo
        enddo
        
        ! Qcr:
        do j=1,NZB+1
           do i=1,NYB+1
              slicedata(i,j,QCRT_TP,slcx_zone)= &
                   real(intFactor_low*Qcr(icellx,i,j)+intFactor_high*Qcr(icellx+1,i,j),KIND=4)
           enddo
        enddo
        
     enddo ! Loop over Blocks

#ifdef ONE_SLC_FILE
     ! Master opens slice file and write its own slice info:
     if (mype .eq. MASTER_PE) then

        slcx_zone_total = 0
        ! Write Block slices into data file:
        do izone = 1,slcx_zone

           slcx_zone_total = slcx_zone_total + 1
           write(index_lb,"(I6.6)") slcx_zone_total
           i = TecZne(                                                &
               'ZONE T=BLKPROC'//index_lb//NULLCHR,                   &
                1,NYB+1,NZB+1,                                        &
               'BLOCK'//NULLCHR,                                      &
                CHAR(0))
 
           ! Write Values:
           i = TecDat(ijkx,slicedata(:,:,XPOS_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,YPOS_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,ZPOS_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,UVEL_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,VVEL_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,WVEL_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,PRES_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,XVOR_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,YVOR_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,ZVOR_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,QCRT_TP,izone),0)
           
        enddo
     endif

   
     ! Now Communicate and write:
     if (mype .eq. MASTER_PE) then

        do iproc = MASTER_PE+1,gr_meshNumProcs-1
           ! Receive number of slices from iproc:
           call mpi_recv(slcx_zone,CONSTANT_ONE,FLASH_INTEGER,iproc,iproc,gr_meshComm,&
                         status,ierr)

           if(slcx_zone .eq. 0) cycle

           ! Receive slice Data:
           buffer_size = (NYB+1)*(NZB+1)*NUM_EXTP*slcx_zone
           call mpi_recv(slicedata(:,:,:,1:slcx_zone),buffer_size,MPI_FLOAT,&
                         iproc,iproc,gr_meshComm,status,ierr)

           ! Write Block Slices into data file:
           do izone = 1 , slcx_zone

              slcx_zone_total = slcx_zone_total + 1
              write(index_lb,"(I6.6)") slcx_zone_total
              i = TecZne(                                                &
                  'ZONE T=BLKPROC'//index_lb//NULLCHR,                   &
                   1,NYB+1,NZB+1,                                        &
                  'BLOCK'//NULLCHR,                                      &
                   CHAR(0))
 
              ! Write Values:
              i = TecDat(ijkx,slicedata(:,:,XPOS_TP,izone),0)
              i = TecDat(ijkx,slicedata(:,:,YPOS_TP,izone),0)
              i = TecDat(ijkx,slicedata(:,:,ZPOS_TP,izone),0)
              i = TecDat(ijkx,slicedata(:,:,UVEL_TP,izone),0)
              i = TecDat(ijkx,slicedata(:,:,VVEL_TP,izone),0)
              i = TecDat(ijkx,slicedata(:,:,WVEL_TP,izone),0)
              i = TecDat(ijkx,slicedata(:,:,PRES_TP,izone),0)
              i = TecDat(ijkx,slicedata(:,:,XVOR_TP,izone),0)
              i = TecDat(ijkx,slicedata(:,:,YVOR_TP,izone),0)
              i = TecDat(ijkx,slicedata(:,:,ZVOR_TP,izone),0)
              i = TecDat(ijkx,slicedata(:,:,QCRT_TP,izone),0)

           enddo
        enddo

     else ! Procs other than MASTER send to it for writing.
  
        ! Send Number of zones:
        call mpi_send(slcx_zone,CONSTANT_ONE,FLASH_INTEGER,MASTER_PE,mype,gr_meshComm,ierr)

        ! Send zones:
        if (slcx_zone .gt. 0) then
           buffer_size = (NYB+1)*(NZB+1)*NUM_EXTP*slcx_zone
           call mpi_send(slicedata(:,:,:,1:slcx_zone),buffer_size,MPI_FLOAT,&
                         MASTER_PE,mype,gr_meshComm,ierr)        
        endif

     endif ! Send and receives


     ! Master closes the file:
     if (mype .eq. MASTER_PE) i = TecEnd()

#else

     if (slcx_zone .gt. 0) then

     write(filename2,'("./IOData/slcx.",i4.4,".",i2.2,".",i6.6,".plt")') &
           count, islicex, mype

     i = TecIni('AMR3D'//NULLCHR,                            &
                'x y z u v w p wx wy wz Q'//NULLCHR,         & !div TV
                filename2//NULLCHR,                          &
                './IOData/'//NULLCHR,                        &
                Debug,VIsdouble)


     slcx_zone_total = 0
     ! Write Block slices into data file:
     do izone = 1,slcx_zone

           slcx_zone_total = slcx_zone_total + 1
           write(index_lb,"(I6.6)") slcx_zone_total
           i = TecZne(                                                &
               'ZONE T=BLKPROC'//index_lb//NULLCHR,                   &
                1,NYB+1,NZB+1,                                        &
               'BLOCK'//NULLCHR,                                      &
                CHAR(0))

           ! Write Values:
           i = TecDat(ijkx,slicedata(:,:,XPOS_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,YPOS_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,ZPOS_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,UVEL_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,VVEL_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,WVEL_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,PRES_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,XVOR_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,YVOR_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,ZVOR_TP,izone),0)
           i = TecDat(ijkx,slicedata(:,:,QCRT_TP,izone),0)

    enddo

    i = TecEnd() 

    endif

#endif


  enddo ! Loop over slices
  deallocate(slicedata)



  ! Write solution to y slices:
  ! write solution data to slcy.XXXX.XX.plt
  allocate(slicedata(NXB+1,NZB+1,NUM_EXTP,MAXBLOCKS)); slicedata = 0.
  do islicey =1,nslicey

     Yslci = Yslice(islicey) 

#ifdef ONE_SLC_FILE
     if (mype .eq. MASTER_PE) then
        write(filename,'("./IOData/slcy.",i4.4,".",i2.2,".plt")') &
              count, islicey

        i = TecIni('AMR3D'//NULLCHR,                            &
                   'x y z u v w p wx wy wz Q'//NULLCHR,         & !div TV
                   filename//NULLCHR,                           &
                   './IOData/'//NULLCHR,                        &
                   Debug,VIsdouble)
     endif
#endif

     ! Loop over blocks:
     slcy_zone = 0
     do lb = 1,blockcount
  
        blockID =  blockList(lb)      

#ifdef FLASH_GRID_PARAMESH
        ! Get Coord and Bsize for the block:
        ! Bounding box:
        call Grid_getBlkBoundBox(blockId,boundBox)
#endif

        ! Test if slice falls on block domain
        if (Yslci .gt. boundBox(2,JAXIS)) cycle
        if (Yslci .le. boundBox(1,JAXIS)) cycle 

#ifdef FLASH_GRID_PARAMESH
        ! Get blocks dx, dy ,dz:
        call Grid_getDeltas(blockID,del)
        dx = del(IAXIS)
        dy = del(JAXIS)
        dz = del(KAXIS)
        
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

        ! V velocity: v(NXB+1,NYB+1,NZB+1)
        ! --------------------------------                           
        tpv = 0.25*(facevaryy(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD:nzc-1) + &
             facevaryy(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD:nzc-1) + &
             facevaryy(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc) + &
             facevaryy(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD+1:nzc));      

        ! W velocity: w(NXB+1,NYB+1,NZB+1)
        ! --------------------------------
        tpw = 0.25*(facevarzz(NGUARD:nxc-1,NGUARD:nyc-1,NGUARD+1:nzc)   + &
             facevarzz(NGUARD+1:nxc,NGUARD:nyc-1,NGUARD+1:nzc)   + &
             facevarzz(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc)   + &
             facevarzz(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD+1:nzc)); 

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
        call ins_velgradtensor(NGUARD,facexData,faceyData,facezData,           &
                               dx,dy,dz,tpdudxc,tpdudyc,tpdudzc,               &
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
        Qcr = 0.5*(2.*(Oxy*Oxy+Oyz*Oyz+Oxz*Oxz) -                &
             (tpdudxcorn*tpdudxcorn + tpdvdycorn*tpdvdycorn +    &
              tpdwdzcorn*tpdwdzcorn + 2.*(Sxy*Sxy+Syz*Syz+Sxz*Sxz)))                      

        ! Release Pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

        ! Interpolate values to slice location:
        icelly = floor((Yslci - yedge(1))/dy) + 1
        if (icelly .eq. (NYB+1)) icelly = icelly - 1

        ! Linear interpolation factors
        intFactor_low  = (yedge(icelly+1)-Yslci)/dy
        intFactor_high = (Yslci - yedge(icelly))/dy

        ! Interpolate to slice array:
        slcy_zone = slcy_zone + 1
        ! X position:
        do j=1,NZB+1
           do i=1,NXB+1
              slicedata(i,j,XPOS_TP,slcy_zone)=real(xedge(i),KIND=4) 
           enddo
        enddo
        
        ! Y Position:
        do j=1,NZB+1
           do i=1,NXB+1
              slicedata(i,j,YPOS_TP,slcy_zone)= &
                   real(intFactor_low*yedge(icelly)+intFactor_high*yedge(icelly+1),KIND=4)
           enddo
        enddo
        
        ! Z Position:
        do j=1,NZB+1
           do i=1,NXB+1
              slicedata(i,j,ZPOS_TP,slcy_zone)=real(zedge(j),KIND=4)
           enddo
        enddo
        
        ! U velocity:
        do j=1,NZB+1
           do i=1,NXB+1
              slicedata(i,j,UVEL_TP,slcy_zone)= &
                   real(intFactor_low*tpu(i,icelly,j)+intFactor_high*tpu(i,icelly+1,j),KIND=4)
           enddo
        enddo
        
        ! V velocity:
        do j=1,NZB+1
           do i=1,NXB+1
              slicedata(i,j,VVEL_TP,slcy_zone)= &
                   real(intFactor_low*tpv(i,icelly,j)+intFactor_high*tpv(i,icelly+1,j),KIND=4)
           enddo
        enddo
        
        ! W velocity:
        do j=1,NZB+1
           do i=1,NXB+1
              slicedata(i,j,WVEL_TP,slcy_zone)= &
                   real(intFactor_low*tpw(i,icelly,j)+intFactor_high*tpw(i,icelly+1,j),KIND=4)
           enddo
        enddo
        
        ! Pressure:
        do j=1,NZB+1
           do i=1,NXB+1
              slicedata(i,j,PRES_TP,slcy_zone)= &
                   real(intFactor_low*tpp(i,icelly,j)+intFactor_high*tpp(i,icelly+1,j),KIND=4)
           enddo
        enddo
        
        ! Vortx:
        do j=1,NZB+1
           do i=1,NXB+1
              slicedata(i,j,XVOR_TP,slcy_zone)= &
                   real(intFactor_low*vortx(i,icelly,j)+intFactor_high*vortx(i,icelly+1,j),KIND=4)
           enddo
        enddo
        
        ! Vorty:
        do j=1,NZB+1
           do i=1,NXB+1
              slicedata(i,j,YVOR_TP,slcy_zone)= &
                   real(intFactor_low*vorty(i,icelly,j)+intFactor_high*vorty(i,icelly+1,j),KIND=4)
           enddo
        enddo
        
        ! Vortz:
        do j=1,NZB+1
           do i=1,NXB+1
              slicedata(i,j,ZVOR_TP,slcy_zone)= &
                   real(intFactor_low*vortz(i,icelly,j)+intFactor_high*vortz(i,icelly+1,j),KIND=4)
           enddo
        enddo
        
        ! Qcr:
        do j=1,NZB+1
           do i=1,NXB+1
              slicedata(i,j,QCRT_TP,slcy_zone)= &
                   real(intFactor_low*Qcr(i,icelly,j)+intFactor_high*Qcr(i,icelly+1,j),KIND=4)
           enddo
        enddo
        
     enddo ! Loop over Blocks

#ifdef ONE_SLC_FILE
     ! Master opens slice file and write its own slice info:
     if (mype .eq. MASTER_PE) then

        slcy_zone_total = 0
        ! Write Block slices into data file:
        do izone = 1,slcy_zone

           slcy_zone_total = slcy_zone_total + 1
           write(index_lb,"(I6.6)") slcy_zone_total
           i = TecZne(                                                &
               'ZONE T=BLKPROC'//index_lb//NULLCHR,                   &
                NXB+1,1,NZB+1,                                        &
               'BLOCK'//NULLCHR,                                      &
                CHAR(0))
 
           ! Write Values:
           i = TecDat(ijky,slicedata(:,:,XPOS_TP,izone),0)
           i = TecDat(ijky,slicedata(:,:,YPOS_TP,izone),0)
           i = TecDat(ijky,slicedata(:,:,ZPOS_TP,izone),0)
           i = TecDat(ijky,slicedata(:,:,UVEL_TP,izone),0)
           i = TecDat(ijky,slicedata(:,:,VVEL_TP,izone),0)
           i = TecDat(ijky,slicedata(:,:,WVEL_TP,izone),0)
           i = TecDat(ijky,slicedata(:,:,PRES_TP,izone),0)
           i = TecDat(ijky,slicedata(:,:,XVOR_TP,izone),0)
           i = TecDat(ijky,slicedata(:,:,YVOR_TP,izone),0)
           i = TecDat(ijky,slicedata(:,:,ZVOR_TP,izone),0)
           i = TecDat(ijky,slicedata(:,:,QCRT_TP,izone),0)
           
        enddo
     endif

   
     ! Now Communicate and write:
     if (mype .eq. MASTER_PE) then

        do iproc = MASTER_PE+1,gr_meshNumProcs-1
           ! Receive number of slices from iproc:
           call mpi_recv(slcy_zone,CONSTANT_ONE,FLASH_INTEGER,iproc,iproc,gr_meshComm,&
                         status,ierr)

           if(slcy_zone .eq. 0) cycle

           ! Receive slice Data:
           buffer_size = (NXB+1)*(NZB+1)*NUM_EXTP*slcy_zone
           call mpi_recv(slicedata(:,:,:,1:slcy_zone),buffer_size,MPI_FLOAT,&
                         iproc,iproc,gr_meshComm,status,ierr)

           ! Write Block Slices into data file:
           do izone = 1 , slcy_zone

              slcy_zone_total = slcy_zone_total + 1
              write(index_lb,"(I6.6)") slcy_zone_total
              i = TecZne(                                                &
                  'ZONE T=BLKPROC'//index_lb//NULLCHR,                   &
                   NXB+1,1,NZB+1,                                        &
                  'BLOCK'//NULLCHR,                                      &
                   CHAR(0))
 
              ! Write Values:
              i = TecDat(ijky,slicedata(:,:,XPOS_TP,izone),0)
              i = TecDat(ijky,slicedata(:,:,YPOS_TP,izone),0)
              i = TecDat(ijky,slicedata(:,:,ZPOS_TP,izone),0)
              i = TecDat(ijky,slicedata(:,:,UVEL_TP,izone),0)
              i = TecDat(ijky,slicedata(:,:,VVEL_TP,izone),0)
              i = TecDat(ijky,slicedata(:,:,WVEL_TP,izone),0)
              i = TecDat(ijky,slicedata(:,:,PRES_TP,izone),0)
              i = TecDat(ijky,slicedata(:,:,XVOR_TP,izone),0)
              i = TecDat(ijky,slicedata(:,:,YVOR_TP,izone),0)
              i = TecDat(ijky,slicedata(:,:,ZVOR_TP,izone),0)
              i = TecDat(ijky,slicedata(:,:,QCRT_TP,izone),0)

           enddo
        enddo

     else ! Procs other than MASTER send to it for writing.
  
        ! Send Number of zones:
        call mpi_send(slcy_zone,CONSTANT_ONE,FLASH_INTEGER,MASTER_PE,mype,gr_meshComm,ierr)

        ! Send zones:
        if (slcy_zone .gt. 0) then
           buffer_size = (NXB+1)*(NZB+1)*NUM_EXTP*slcy_zone
           call mpi_send(slicedata(:,:,:,1:slcy_zone),buffer_size,MPI_FLOAT,&
                         MASTER_PE,mype,gr_meshComm,ierr)        
        endif

     endif ! Send and receives


     ! Master closes the file:
     if (mype .eq. MASTER_PE) i = TecEnd()

#else

     if (slcy_zone .gt. 0) then

     write(filename2,'("./IOData/slcy.",i4.4,".",i2.2,".",i6.6,".plt")') &
           count, islicey, mype

     i = TecIni('AMR3D'//NULLCHR,                            &
                'x y z u v w p wx wy wz Q'//NULLCHR,         & !div TV
                filename2//NULLCHR,                          &
                './IOData/'//NULLCHR,                        &
                Debug,VIsdouble)

     slcy_zone_total = 0
     ! Write Block slices into data file:
     do izone = 1,slcy_zone

        slcy_zone_total = slcy_zone_total + 1
        write(index_lb,"(I6.6)") slcy_zone_total
        i = TecZne(                                                &
            'ZONE T=BLKPROC'//index_lb//NULLCHR,                   &
             NXB+1,1,NZB+1,                                        &
            'BLOCK'//NULLCHR,                                      &
             CHAR(0))

        ! Write Values:
        i = TecDat(ijky,slicedata(:,:,XPOS_TP,izone),0)
        i = TecDat(ijky,slicedata(:,:,YPOS_TP,izone),0)
        i = TecDat(ijky,slicedata(:,:,ZPOS_TP,izone),0)
        i = TecDat(ijky,slicedata(:,:,UVEL_TP,izone),0)
        i = TecDat(ijky,slicedata(:,:,VVEL_TP,izone),0)
        i = TecDat(ijky,slicedata(:,:,WVEL_TP,izone),0)
        i = TecDat(ijky,slicedata(:,:,PRES_TP,izone),0)
        i = TecDat(ijky,slicedata(:,:,XVOR_TP,izone),0)
        i = TecDat(ijky,slicedata(:,:,YVOR_TP,izone),0)
        i = TecDat(ijky,slicedata(:,:,ZVOR_TP,izone),0)
        i = TecDat(ijky,slicedata(:,:,QCRT_TP,izone),0)

     enddo

     i = TecEnd()

     endif

#endif

  enddo ! Loop over slices
  deallocate(slicedata)

  ! Write solution to z slices:
  ! write solution data to slcz.XXXX.XX.plt
  allocate(slicedata(NXB+1,NYB+1,NUM_EXTP,MAXBLOCKS)); slicedata = 0.
  do islicez =1,nslicez

     Zslci = Zslice(islicez) 

#ifdef ONE_SLC_FILE
     if (mype .eq. MASTER_PE) then
        write(filename,'("./IOData/slcz.",i4.4,".",i2.2,".plt")') &
              count, islicez

        i = TecIni('AMR3D'//NULLCHR,                            &
                   'x y z u v w p wx wy wz Q'//NULLCHR,         & !div TV
                   filename//NULLCHR,                           &
                   './IOData/'//NULLCHR,                        &
                   Debug,VIsdouble)
     endif
#endif

     ! Loop over blocks:
     slcz_zone = 0
     do lb = 1,blockcount
  
        blockID =  blockList(lb)      

#ifdef FLASH_GRID_PARAMESH
        ! Get Coord and Bsize for the block:
        ! Bounding box:
        call Grid_getBlkBoundBox(blockId,boundBox)
#endif

        ! Test if slice falls on block domain
        if (Zslci .gt. boundBox(2,KAXIS)) cycle
        if (Zslci .le. boundBox(1,KAXIS)) cycle 

#ifdef FLASH_GRID_PARAMESH
        ! Get blocks dx, dy ,dz:
        call Grid_getDeltas(blockID,del)
        dx = del(IAXIS)
        dy = del(JAXIS)
        dz = del(KAXIS)
        
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

        ! V velocity: v(NXB+1,NYB+1,NZB+1)
        ! --------------------------------                           
        tpv = 0.25*(facevaryy(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD:nzc-1) + &
             facevaryy(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD:nzc-1) + &
             facevaryy(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc) + &
             facevaryy(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD+1:nzc));      

        ! W velocity: w(NXB+1,NYB+1,NZB+1)
        ! --------------------------------
        tpw = 0.25*(facevarzz(NGUARD:nxc-1,NGUARD:nyc-1,NGUARD+1:nzc)   + &
             facevarzz(NGUARD+1:nxc,NGUARD:nyc-1,NGUARD+1:nzc)   + &
             facevarzz(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc)   + &
             facevarzz(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD+1:nzc)); 

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
        call ins_velgradtensor(NGUARD,facexData,faceyData,facezData,           &
                               dx,dy,dz,tpdudxc,tpdudyc,tpdudzc,               &
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
        Qcr = 0.5*(2.*(Oxy*Oxy+Oyz*Oyz+Oxz*Oxz) -                &
             (tpdudxcorn*tpdudxcorn + tpdvdycorn*tpdvdycorn +    &
              tpdwdzcorn*tpdwdzcorn + 2.*(Sxy*Sxy+Syz*Syz+Sxz*Sxz)))                      

        ! Release Pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

        ! Interpolate values to slice location:
        icellz = floor((Zslci - zedge(1))/dz) + 1
        if (icellz .eq. (NZB+1)) icellz = icellz - 1

        ! Linear interpolation factors
        intFactor_low  = (zedge(icellz+1)-Zslci)/dz
        intFactor_high = (Zslci - zedge(icellz))/dz

        ! Interpolate to slice array:
        slcz_zone = slcz_zone + 1
        ! X position:
        do j=1,NYB+1
           do i=1,NXB+1
              slicedata(i,j,XPOS_TP,slcz_zone)=real(xedge(i),KIND=4)
           enddo
        enddo
        
        ! Y Position:
        do j=1,NYB+1
           do i=1,NXB+1
              slicedata(i,j,YPOS_TP,slcz_zone)=real(yedge(j),KIND=4)
           enddo
        enddo
        
        ! Z Position:
        do j=1,NYB+1
           do i=1,NXB+1
              slicedata(i,j,ZPOS_TP,slcz_zone)= &
                   real(intFactor_low*zedge(icellz)+intFactor_high*zedge(icellz+1),KIND=4)
           enddo
        enddo
        
        ! U velocity:
        do j=1,NYB+1
           do i=1,NXB+1
              slicedata(i,j,UVEL_TP,slcz_zone)= &
                   real(intFactor_low*tpu(i,j,icellz)+intFactor_high*tpu(i,j,icellz+1),KIND=4)
           enddo
        enddo
        
        ! V velocity:
        do j=1,NYB+1
           do i=1,NXB+1
              slicedata(i,j,VVEL_TP,slcz_zone)= &
                   real(intFactor_low*tpv(i,j,icellz)+intFactor_high*tpv(i,j,icellz+1),KIND=4)
           enddo
        enddo
        
        ! W velocity:
        do j=1,NYB+1
           do i=1,NXB+1
              slicedata(i,j,WVEL_TP,slcz_zone)= &
                   real(intFactor_low*tpw(i,j,icellz)+intFactor_high*tpw(i,j,icellz+1),KIND=4)
           enddo
        enddo
        
        ! Pressure:
        do j=1,NYB+1
           do i=1,NXB+1
              slicedata(i,j,PRES_TP,slcz_zone)= &
                   real(intFactor_low*tpp(i,j,icellz)+intFactor_high*tpp(i,j,icellz+1),KIND=4)
           enddo
        enddo
        
        ! Vortx:
        do j=1,NYB+1
           do i=1,NXB+1
              slicedata(i,j,XVOR_TP,slcz_zone)= &
                   real(intFactor_low*vortx(i,j,icellz)+intFactor_high*vortx(i,j,icellz+1),KIND=4)
           enddo
        enddo
        
        ! Vorty:
        do j=1,NYB+1
           do i=1,NXB+1
              slicedata(i,j,YVOR_TP,slcz_zone)= &
                   real(intFactor_low*vorty(i,j,icellz)+intFactor_high*vorty(i,j,icellz+1),KIND=4)
           enddo
        enddo
        
        ! Vortz:
        do j=1,NYB+1
           do i=1,NXB+1
              slicedata(i,j,ZVOR_TP,slcz_zone)= &
                   real(intFactor_low*vortz(i,j,icellz)+intFactor_high*vortz(i,j,icellz+1),KIND=4)
           enddo
        enddo
        
        ! Qcr:
        do j=1,NYB+1
           do i=1,NXB+1
              slicedata(i,j,QCRT_TP,slcz_zone)= &
                   real(intFactor_low*Qcr(i,j,icellz)+intFactor_high*Qcr(i,j,icellz+1),KIND=4)
           enddo
        enddo
        
     enddo ! Loop over Blocks


#ifdef ONE_SLC_FILE

     ! Master opens slice file and write its own slice info:
     if (mype .eq. MASTER_PE) then

        slcz_zone_total = 0
        ! Write Block slices into data file:
        do izone = 1,slcz_zone

           slcz_zone_total = slcz_zone_total + 1
           write(index_lb,"(I6.6)") slcz_zone_total
           i = TecZne(                                                &
               'ZONE T=BLKPROC'//index_lb//NULLCHR,                   &
                NXB+1,NYB+1,1,                                        &
               'BLOCK'//NULLCHR,                                      &
                CHAR(0))
 
           ! Write Values:
           i = TecDat(ijkz,slicedata(:,:,XPOS_TP,izone),0)
           i = TecDat(ijkz,slicedata(:,:,YPOS_TP,izone),0)
           i = TecDat(ijkz,slicedata(:,:,ZPOS_TP,izone),0)
           i = TecDat(ijkz,slicedata(:,:,UVEL_TP,izone),0)
           i = TecDat(ijkz,slicedata(:,:,VVEL_TP,izone),0)
           i = TecDat(ijkz,slicedata(:,:,WVEL_TP,izone),0)
           i = TecDat(ijkz,slicedata(:,:,PRES_TP,izone),0)
           i = TecDat(ijkz,slicedata(:,:,XVOR_TP,izone),0)
           i = TecDat(ijkz,slicedata(:,:,YVOR_TP,izone),0)
           i = TecDat(ijkz,slicedata(:,:,ZVOR_TP,izone),0)
           i = TecDat(ijkz,slicedata(:,:,QCRT_TP,izone),0)
           
        enddo
     endif

   
     ! Now Communicate and write:
     if (mype .eq. MASTER_PE) then

        do iproc = MASTER_PE+1,gr_meshNumProcs-1
           ! Receive number of slices from iproc:
           call mpi_recv(slcz_zone,CONSTANT_ONE,FLASH_INTEGER,iproc,iproc,gr_meshComm,&
                         status,ierr)

           if(slcz_zone .eq. 0) cycle

           ! Receive slice Data:
           buffer_size = (NXB+1)*(NYB+1)*NUM_EXTP*slcz_zone
           call mpi_recv(slicedata(:,:,:,1:slcz_zone),buffer_size,MPI_FLOAT,&
                         iproc,iproc,gr_meshComm,status,ierr)

           ! Write Block Slices into data file:
           do izone = 1 , slcz_zone

              slcz_zone_total = slcz_zone_total + 1
              write(index_lb,"(I6.6)") slcz_zone_total
              i = TecZne(                                                &
                  'ZONE T=BLKPROC'//index_lb//NULLCHR,                   &
                   NXB+1,NYB+1,1,                                        &
                  'BLOCK'//NULLCHR,                                      &
                   CHAR(0))
 
              ! Write Values:
              i = TecDat(ijkz,slicedata(:,:,XPOS_TP,izone),0)
              i = TecDat(ijkz,slicedata(:,:,YPOS_TP,izone),0)
              i = TecDat(ijkz,slicedata(:,:,ZPOS_TP,izone),0)
              i = TecDat(ijkz,slicedata(:,:,UVEL_TP,izone),0)
              i = TecDat(ijkz,slicedata(:,:,VVEL_TP,izone),0)
              i = TecDat(ijkz,slicedata(:,:,WVEL_TP,izone),0)
              i = TecDat(ijkz,slicedata(:,:,PRES_TP,izone),0)
              i = TecDat(ijkz,slicedata(:,:,XVOR_TP,izone),0)
              i = TecDat(ijkz,slicedata(:,:,YVOR_TP,izone),0)
              i = TecDat(ijkz,slicedata(:,:,ZVOR_TP,izone),0)
              i = TecDat(ijkz,slicedata(:,:,QCRT_TP,izone),0)

           enddo
        enddo

     else ! Procs other than MASTER send to it for writing.
  
        ! Send Number of zones:
        call mpi_send(slcz_zone,CONSTANT_ONE,FLASH_INTEGER,MASTER_PE,mype,gr_meshComm,ierr)

        ! Send zones:
        if (slcz_zone .gt. 0) then
           buffer_size = (NXB+1)*(NYB+1)*NUM_EXTP*slcz_zone
           call mpi_send(slicedata(:,:,:,1:slcz_zone),buffer_size,MPI_FLOAT,&
                         MASTER_PE,mype,gr_meshComm,ierr)        
        endif

     endif ! Send and receives


     ! Master closes the file:
     if (mype .eq. MASTER_PE) i = TecEnd()

#else


     if (slcz_zone .gt. 0) then

     write(filename2,'("./IOData/slcz.",i4.4,".",i2.2,".",i6.6,".plt")') &
           count, islicez, mype

     i = TecIni('AMR3D'//NULLCHR,                            &
                'x y z u v w p wx wy wz Q'//NULLCHR,         & !div TV
                filename2//NULLCHR,                          &
                './IOData/'//NULLCHR,                        &
                Debug,VIsdouble)

     slcz_zone_total = 0
     ! Write Block slices into data file:
     do izone = 1,slcz_zone

        slcz_zone_total = slcz_zone_total + 1
        write(index_lb,"(I6.6)") slcz_zone_total
        i = TecZne(                                                &
            'ZONE T=BLKPROC'//index_lb//NULLCHR,                   &
             NXB+1,NYB+1,1,                                        &
            'BLOCK'//NULLCHR,                                      &
             CHAR(0))

        ! Write Values:
        i = TecDat(ijkz,slicedata(:,:,XPOS_TP,izone),0)
        i = TecDat(ijkz,slicedata(:,:,YPOS_TP,izone),0)
        i = TecDat(ijkz,slicedata(:,:,ZPOS_TP,izone),0)
        i = TecDat(ijkz,slicedata(:,:,UVEL_TP,izone),0)
        i = TecDat(ijkz,slicedata(:,:,VVEL_TP,izone),0)
        i = TecDat(ijkz,slicedata(:,:,WVEL_TP,izone),0)
        i = TecDat(ijkz,slicedata(:,:,PRES_TP,izone),0)
        i = TecDat(ijkz,slicedata(:,:,XVOR_TP,izone),0)
        i = TecDat(ijkz,slicedata(:,:,YVOR_TP,izone),0)
        i = TecDat(ijkz,slicedata(:,:,ZVOR_TP,izone),0)
        i = TecDat(ijkz,slicedata(:,:,QCRT_TP,izone),0)

     enddo

     i = TecEnd()

     endif


#endif

  enddo ! Loop over slices
  deallocate(slicedata)

#ifdef ALLOC_TECPLOT_VARS

  deallocate(facevarxx,facevaryy,facevarzz)
  deallocate(tpu,tpv,tpw,tpp)
  deallocate(tpdudxcorn, tpdudycorn, tpdudzcorn)
  deallocate(tpdvdxcorn, tpdvdycorn, tpdvdzcorn)
  deallocate(tpdwdxcorn, tpdwdycorn, tpdwdzcorn)
  deallocate(vortx,vorty,vortz,omg)
  deallocate(Sxy,Syz,Sxz,Oxy,Oyz,Oxz,Qcr,divpp,TVtpp)
  deallocate(arraylb)
  deallocate(tpdudxc,tpdudyc,tpdudzc,tpdvdxc,tpdvdyc,tpdvdzc,tpdwdxc,tpdwdyc,tpdwdzc)

#endif

  if (mype .eq. 0) then

    if (nslicex .gt. 0) then
     write(*,*) ''
     write(filename,'("./IOData/slcx.",i4.4,".**.plt")') count
     write(*,*) '*** Wrote plotfile to ',filename,' ****'
    endif

    if (nslicey .gt. 0) then
     write(filename,'("./IOData/slcy.",i4.4,".**.plt")') count
     write(*,*) '*** Wrote plotfile to ',filename,' ****'
    endif

    if (nslicez .gt. 0) then
     write(filename,'("./IOData/slcz.",i4.4,".**.plt")') count
     write(*,*) '*** Wrote plotfile to ',filename,' ****'
    endif

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

