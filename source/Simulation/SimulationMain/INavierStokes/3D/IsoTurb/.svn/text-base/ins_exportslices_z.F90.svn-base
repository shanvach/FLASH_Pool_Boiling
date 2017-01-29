! subroutine ins_exportslices_z
! Explorts slices on Z planes for turbulence simulations:
! One file will be written per slice.
! Values are exported at cell centers.
! --------------------------------------------------------------------------

  subroutine ins_exportslices_z(istep,iwrite,time,mvisc,blockCount,blockList)

  ! Modules Use:
  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkCenterCoords,&
                             Grid_getBlkPhysicalSize,&
                             Grid_getBlkBoundBox,    &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_getBlkCornerID,    &
                             Grid_getGlobalIndexLimits

  
  use ins_interface, only  :  ins_velomg2center
  
  use Grid_data, only : gr_meshMe, gr_meshNumProcs, gr_meshComm, gr_imin, gr_imax, gr_jmin, gr_jmax
  
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "IncompNS.h"

  integer :: istep, iwrite
  real    :: time, mvisc
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(blockCount) :: blockList

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
  
  real :: del(MDIM),coord(MDIM),bsize(MDIM),boundBox(CONSTANT_TWO,MDIM),dx,dy,dz
  real :: volcell,tvol(MAXBLOCKS)

  ! Variables on slices, function of z:

!!$#ifdef COMPUTE_STATS
!!$
!!$  ! Index parameters:
!!$  integer, parameter :: UAVG_IND  = 1
!!$  integer, parameter :: VAVG_IND  = 2
!!$  integer, parameter :: WAVG_IND  = 3
!!$  integer, parameter :: OMGX2_IND = 4
!!$  integer, parameter :: OMGY2_IND = 5
!!$  integer, parameter :: OMGZ2_IND = 6 
!!$  integer, parameter :: UU_IND    = 7
!!$  integer, parameter :: VV_IND    = 8
!!$  integer, parameter :: WW_IND    = 9
!!$
!!$  ! Plane Average squared vorticites vs. z:
!!$  real :: omgx2(NZB,MAXBLOCKS),omgy2(NZB,MAXBLOCKS),omgz2(NZB,MAXBLOCKS)
!!$  ! Plane average Reynolds stresses vs. z:
!!$  real :: uu(NZB,MAXBLOCKS),vv(NZB,MAXBLOCKS),ww(NZB,MAXBLOCKS), &
!!$          uv(NZB,MAXBLOCKS),uw(NZB,MAXBLOCKS),vw(NZB,MAXBLOCKS)
!!$  ! Plane average q^2 -> TKE = 1/2 * q^2  vs. z:
!!$  real :: q2(NZB,MAXBLOCKS)
!!$  ! Plane average of rms fluctuations vs. z:
!!$  real :: urms2(NZB,MAXBLOCKS),vrms2(NZB,MAXBLOCKS),wrms2(NZB,MAXBLOCKS)
!!$  ! Plane average of longitudinal and transverse derivatives vs. z:
!!$  real :: dudx2(NZB,MAXBLOCKS),dudy2(NZB,MAXBLOCKS)
!!$  ! Plane average of longitudinal and transverse microscales vs. z:
!!$  real :: lambdaf(NZB,MAXBLOCKS),lambdag(NZB,MAXBLOCKS)
!!$  ! Plane average of dissipation vs. z:
!!$  real :: eps_nu(NZB,MAXBLOCKS)
!!$
!!$  ! If LES
!!$  ! Plane average of turbulent viscosity and SGS dissipation vs. z:
!!$  real :: nut(NZB,MAXBLOCKS),eps_nut(NZB,MAXBLOCKS)
!!$
!!$  ! Average presssure field, and rms fluctuation:
!!$  real :: pavg(NZB,MAXBLOCKS),prms2(NZB,MAXBLOCKS)
!!$
!!$  real :: dux,dvy,dwz,duy,duz,dvx,dvz,dwx,dwy
!!$
!!$
!!$#else /* write out variables only */

  ! Instantaneous field indexes
  integer, parameter :: U_IND     = 1
  integer, parameter :: V_IND     = 2
  integer, parameter :: W_IND     = 3
  
  integer, parameter :: DUDX_IND  = 4    
  integer, parameter :: DUDY_IND  = 5   
  integer, parameter :: DUDZ_IND  = 6   
  integer, parameter :: DVDX_IND  = 7    
  integer, parameter :: DVDY_IND  = 8    
  integer, parameter :: DVDZ_IND  = 9    
  integer, parameter :: DWDX_IND  =10    
  integer, parameter :: DWDY_IND  =11    
  integer, parameter :: DWDZ_IND  =12   

  integer, parameter :: PRES_IND  =13
  integer, parameter :: TVIS_IND  =14

  ! Here indexes for multiphase rhos, phi etc.


  ! Total vars:
  integer, parameter :: TotVars   =14

!!$#endif

  integer, save :: nslicez
  real, save, allocatable, dimension(:) :: zc, Zslice
  logical, save :: readfilez
  real :: st_zmin,st_zmax

  ! Maximum number of points in x and y for a slice:
  integer, parameter :: NXBMAX = 386
  integer, parameter :: NYBMAX = 386

  ! Average velocities on each plane zc:
  real, allocatable, dimension(:,:,:) :: data,dataaux

  integer i,j,k,kk,lb,blockID,ierr

  logical, save :: initfile  = .TRUE.
  logical, save :: firstcall = .TRUE.

  character(80) :: filename

  integer, dimension(MDIM) :: globalIndexLimits,cornerID,stride,strideproc

  real xedge(NXB+1),xcell(NXB+1)
  real yedge(NYB+1),ycell(NYB+1)
  real zedge(NZB+1),zcell(NZB+1)
  real intsx(NXB+1), intsy(NYB+1), intsz(NZB+1)

  integer :: icellz,istartx, istarty, ii, jj, slceproc, totslceproc, islicez

  real :: Zslci, totzc
  logical :: izc

  real, parameter :: eps = 1.e-10 

  !--- Initialize averages:
  allocate(data(NXBMAX,NYBMAX,TotVars))
  allocate(dataaux(NXBMAX,NYBMAX,TotVars))

  ! Call filling of cell centered velocities and vorticity:
  call ins_velomg2center(blockList,blockCount)


  ! Load slice positions
  if (firstcall) then

     call RuntimeParameters_get("zmin", st_zmin)
     call RuntimeParameters_get("zmax", st_zmax)

    ! Z slices:
     call RuntimeParameters_get("stats_nslicesz", nslicez)

     !write(*,*) 'zmin max, nslicez=',st_zmin,st_zmin,nslicez

     if (nslicez .gt. 0) then
        allocate(Zslice(nslicez))
        allocate(zc(nslicez))
        call RuntimeParameters_get("stats_slicesz_from_file", readfilez)
        if (readfilez) then
           open(unit=33, file='./zlist_stats.slc', status='old')
           do islicez = 1,nslicez
              read(33,*) Zslice(islicez)
           end do
           close(33)
        else
           ! write evenly spaced across the domain:
           dz = (st_zmax-st_zmin)/real(nslicez)
           do islicez =1,nslicez
              Zslice(islicez) = st_zmin + 0.5*dz + real(islicez-1)*dz
           enddo
        endif
     endif

     firstcall = .false.

     !write(*,*) 'Zslice=',Zslice

  endif


  ! Now all processors compute and interpolate to the slice:
  intsx    = (/ (real(i), i=0,NXB) /)
  intsy    = (/ (real(i), i=0,NYB) /)
  intsz    = (/ (real(i), i=0,NZB) /)


  ! Loop over slices:
  do islicez = 1 , nslicez

     Zslci = Zslice(islicez)

     zc(islicez)    = 0.
     data = 0.
     izc  = .true.
     slceproc       = 0    
     stride(1:MDIM) = 0

     ! Loop over blocks:
     do lb = 1,blockCount

        blockID = blockList(lb)     

        ! Bounding box:
        call Grid_getBlkBoundBox(blockId,boundBox)

        ! Test if slice falls on block domain
        if (Zslci .gt. boundBox(2,KAXIS)) cycle
        if (Zslci .le. boundBox(1,KAXIS)) cycle


        ! Get blocks dx, dy ,dz, coord, bsize:
        call Grid_getDeltas(blockID,del)
        dx = del(IAXIS)
        dy = del(JAXIS)
        dz = del(KAXIS)

        call Grid_getBlkCenterCoords(blockID,coord)

        bsize(:) = boundBox(HIGH,:) - boundBox(LOW,:)

        ! Cell, face locations:
        xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
        xcell = xedge(:) + dx/2.0;

        yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
        ycell = yedge(:) + dy/2.0;

        zedge = coord(KAXIS) - bsize(KAXIS)/2.0 + dz*intsz;
        zcell = zedge(:) + dz/2.0;


        ! Interpolate values to slice location:
        icellz = floor((Zslci - zedge(1))/dz) + 1
        if (icellz .eq. (NZB+1)) icellz = icellz - 1

        ! Linear interpolation factors
        !intFactor_low  = (zedge(icellz+1)-Zslci)/dz
        !intFactor_high = (Zslci - zedge(icellz))/dz
        
        ! Cell location in z:
        if (izc) then 
           zc(islicez) = zcell(icellz)         
           izc         = .false.
           slceproc    = 1
        else
           if (abs(zc(islicez)-zcell(icellz)) .gt. eps) then

              write(*,*) ' '
              write(*,*) 'MyPe=',gr_meshMe,', block',blockID,coord(:)
              write(*,*) 'zc block=',zcell(icellz),', zc prev block=',zc(islicez)
              write(*,*) 'Slice location=',Zslci 

              call Driver_abortFlash("Z location for blocks on same slice is not the same.")

           endif
        endif

        ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
        call Grid_getBlkPtr(blockID,facezData,FACEZ)   

        ! Block Corner ID
        call Grid_getBlkCornerID(blockId, cornerID, stride)

        ! Indexes for data in x, y:        
        istartx = (cornerID(IAXIS)-1)/stride(IAXIS)
        istarty = (cornerID(JAXIS)-1)/stride(JAXIS)

        ! These are the index limits at maximum refinement: 
        call Grid_getGlobalIndexLimits(globalIndexLimits)
        if ( (istartx+NXB) .gt. NXBMAX) then
             write(*,*) ' '
             write(*,*) 'MyPe=',gr_meshMe,', block',blockID,', max ind X=',(istartx+NXB),NXBMAX
             call Driver_abortFlash("Max index required for data greater than NXBMAX. Increase NXBMAX.")
        endif
        if ( (istarty+NYB) .gt. NYBMAX) then
             write(*,*) ' '
             write(*,*) 'MyPe=',gr_meshMe,', block',blockID,', max ind Y=',(istarty+NYB),NYBMAX
             call Driver_abortFlash("Max index required for data greater than NYBMAX. Increase NYBMAX.")
        endif

        ! Compute variables, load to data:
        jj = istarty
        k  = icellz + NGUARD
        do j = GRID_JLO , GRID_JHI
          jj = jj + 1
          ii = istartx
          do i = GRID_ILO , GRID_IHI
            ii = ii + 1
    
            ! Velocities: From the block face in the collocated direction to
            ! the last internal face in that direction, that is shifted 1/2*del.
            !data(ii,jj,U_IND) = facexData(VELC_FACE_VAR,i,j,k)

            !data(ii,jj,V_IND) = faceyData(VELC_FACE_VAR,i,j,k)

            !data(ii,jj,W_IND) = facezData(VELC_FACE_VAR,i,j,k)

            ! Acerage to cell centers:
            data(ii,jj,U_IND) = 0.5*(facexData(VELC_FACE_VAR,i,j,k)+facexData(VELC_FACE_VAR,i+1,j,k))

            data(ii,jj,V_IND) = 0.5*(faceyData(VELC_FACE_VAR,i,j,k)+faceyData(VELC_FACE_VAR,i,j+1,k))

            data(ii,jj,W_IND) = 0.5*(facezData(VELC_FACE_VAR,i,j,k)+facezData(VELC_FACE_VAR,i,j,k+1))
            


            ! Velocity derivatives:
            data(ii,jj,DUDX_IND)=(facexData(VELC_FACE_VAR,i+1,j,k)-facexData(VELC_FACE_VAR,i,j,k))/del(IAXIS)
            data(ii,jj,DVDY_IND)=(faceyData(VELC_FACE_VAR,i,j+1,k)-faceyData(VELC_FACE_VAR,i,j,k))/del(JAXIS)
            data(ii,jj,DWDZ_IND)=(facezData(VELC_FACE_VAR,i,j,k+1)-facezData(VELC_FACE_VAR,i,j,k))/del(KAXIS)

            data(ii,jj,DUDY_IND)=0.25*(facexData(VELC_FACE_VAR,i+1,j+1,k)-facexData(VELC_FACE_VAR,i+1,j-1,k) + &
                                       facexData(VELC_FACE_VAR,i,j+1,k)  -facexData(VELC_FACE_VAR,i,j-1,k) )/del(JAXIS)
            data(ii,jj,DUDZ_IND)=0.25*(facexData(VELC_FACE_VAR,i+1,j,k+1)-facexData(VELC_FACE_VAR,i+1,j,k-1) + &
                                       facexData(VELC_FACE_VAR,i,j,k+1)  -facexData(VELC_FACE_VAR,i,j,k-1) )/del(KAXIS)

            data(ii,jj,DVDX_IND)=0.25*(faceyData(VELC_FACE_VAR,i+1,j+1,k)-faceyData(VELC_FACE_VAR,i-1,j+1,k) + &  
                                       faceyData(VELC_FACE_VAR,i+1,j,k)  -faceyData(VELC_FACE_VAR,i-1,j,k) )/del(IAXIS)
            data(ii,jj,DVDZ_IND)=0.25*(faceyData(VELC_FACE_VAR,i,j+1,k+1)-faceyData(VELC_FACE_VAR,i,j+1,k-1) + &
                                       faceyData(VELC_FACE_VAR,i,j,k+1)  -faceyData(VELC_FACE_VAR,i,j,k-1) )/del(KAXIS)
            
            data(ii,jj,DWDX_IND)=0.25*(facezData(VELC_FACE_VAR,i+1,j,k+1)-facezData(VELC_FACE_VAR,i-1,j,k+1) + &
                                       facezData(VELC_FACE_VAR,i+1,j,k)  -facezData(VELC_FACE_VAR,i-1,j,k) )/del(IAXIS)
            data(ii,jj,DWDY_IND)=0.25*(facezData(VELC_FACE_VAR,i,j+1,k+1)-facezData(VELC_FACE_VAR,i,j-1,k+1) + &
                                       facezData(VELC_FACE_VAR,i,j+1,k)  -facezData(VELC_FACE_VAR,i,j-1,k) )/del(JAXIS)

            data(ii,jj,PRES_IND)=solnData(PRES_VAR,i,j,k)

            data(ii,jj,TVIS_IND)=solnData(TVIS_VAR,i,j,k)

     
          enddo
        enddo

     enddo ! loop blocks

  
     ! Allreduce zc location and test:     
     ! First read number of processor that contain slice, and their zc(value)
     call mpi_reduce(slceproc, totslceproc, 1, FLASH_INTEGER, MPI_SUM, MASTER_PE, gr_meshComm, ierr)
     call mpi_reduce(zc(islicez), totzc, 1, FLASH_REAL, MPI_SUM, MASTER_PE, gr_meshComm, ierr)  

     ! Now master tests location:
     if (gr_meshME .eq. MASTER_PE) then

        if (totslceproc .eq. 0) then
           write(*,*) ' '
           write(*,*) 'Tot number of procs that contain slice ',Zslci,'is zero.'
     
           call Driver_abortFlash("No processors contain slice location.")

        endif
       
        totzc = totzc/real(totslceproc)

        if ( real(slceproc)*abs(zc(islicez)-totzc) .gt. eps) then ! One of the slice locations doesn't match, abort:

           write(*,*) ' '
           write(*,*) 'MyPe=',gr_meshMe,', slice=',islicez
           write(*,*) 'My zc=',zc(islicez),', doesnt match totzc=',totzc,slceproc

           call Driver_abortFlash("Z location for slice is not the same on different processors.")

        endif

     endif

     ! Compute the NXBtot, NYBtot:
     ! Reduce max stride:
     strideproc(1:MDIM) = stride(1:MDIM)
     call mpi_reduce(strideproc, stride, MDIM, FLASH_INTEGER, MPI_MAX, MASTER_PE, gr_meshComm, ierr)

     if (gr_meshME .eq. MASTER_PE) then

        if (any(stride(1:MDIM) .eq. 0)) then

           write(*,*) ' ' 
           write(*,*) 'Zero stride=',stride(1:MDIM)

           call Driver_abortFlash("Zero stride for slice, slice may be out of domain.")

        endif 
     
        ! These are the index limits at maximum refinement: 
        call Grid_getGlobalIndexLimits(globalIndexLimits)
        ! These are the index limits at the slice refinement level (can be less than max,
        ! but the same level has to be present for the whole slice):
        globalIndexLimits(IAXIS) = globalIndexLimits(IAXIS)/stride(IAXIS)
        globalIndexLimits(JAXIS) = globalIndexLimits(JAXIS)/stride(JAXIS)

     endif

     ! Now Brute force reduce data:
     dataaux = data
     call mpi_reduce(dataaux, data, NXBMAX*NYBMAX*TotVars, FLASH_REAL, MPI_SUM, MASTER_PE, gr_meshComm, ierr)


     ! Master Processor Writes:
     if (gr_meshME .eq. MASTER_PE) then

        ! Slice name:
        write(filename,'("./IOData/stat_slcz.",i4.4,".",i3.3,".bin")') &
                        iwrite, islicez     

        write(*,*) 'Writing filename=',trim(filename),totzc

        open(1,file=trim(filename),status='replace',form='unformatted')

        write(1) time,mvisc,totzc,Zslci
        write(1) istep,globalIndexLimits(IAXIS),globalIndexLimits(JAXIS),TotVars
        write(1) gr_imin,gr_imax,gr_jmin,gr_jmax

        do k = 1, Totvars 

           write(1) ( (data(i,j,k),i=1,globalIndexLimits(IAXIS)) , j=1,globalIndexLimits(JAXIS) )

        enddo

        close(1)

     endif


  enddo ! End do slices:


  deallocate(data,dataaux)

  return
  end subroutine



