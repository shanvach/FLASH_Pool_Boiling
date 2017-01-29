! subroutine ins_turbstats_z
! Computes turbulent statistics on Z planes for turbulence simulations:
! One file will be written per processor, the plane averages done on the
! Processors portion of the plane.
! This routine is intended for the domain in x and y to be completely owned
! by the processor
! Values are exported at cell centers.
! --------------------------------------------------------------------------

  subroutine ins_turbstats_z(istep,iwrite,time,mvisc,blockCount,blockList)

  ! Modules Use:
  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkCenterCoords,&
                             Grid_getBlkPhysicalSize,&
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits

  
  use ins_interface, only  :  ins_velomg2center
  
  use Grid_data, only : gr_meshMe, gr_meshNumProcs
  
  implicit none
#include "constants.h"
#include "Flash.h"
#include "IncompNS.h"
  include 'mpif.h'

  integer :: istep, iwrite
  real    :: time, mvisc
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(blockCount) :: blockList

!  integer, parameter ::  ng = NGUARD
!  integer, parameter ::  nxi= NGUARD + NXB
!  integer, parameter ::  nyj= NGUARD + NYB
!  integer, parameter ::  nzk= NGUARD + NZB
!  integer, parameter ::  nxc= NGUARD + NXB + 1
!  integer, parameter ::  nyc= NGUARD + NYB + 1
!  integer, parameter ::  nzc= NGUARD + NZB + 1

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
  
  real :: del(MDIM),coord(MDIM),bsize(MDIM),dz !,dx,dy,dz
  real :: volcell,tvol(MAXBLOCKS)

  ! Variables to plane average and export, function of z:
  ! Slice z Location:
  real :: zc(NZB,MAXBLOCKS)
  ! Average velocities on each plane vs. z:
  real :: uavg(NZB,MAXBLOCKS),vavg(NZB,MAXBLOCKS),wavg(NZB,MAXBLOCKS)
  ! Plane Average squared vorticites vs. z:
  real :: omgx2(NZB,MAXBLOCKS),omgy2(NZB,MAXBLOCKS),omgz2(NZB,MAXBLOCKS)
  ! Plane average Reynolds stresses vs. z:
  real :: uu(NZB,MAXBLOCKS),vv(NZB,MAXBLOCKS),ww(NZB,MAXBLOCKS), &
          uv(NZB,MAXBLOCKS),uw(NZB,MAXBLOCKS),vw(NZB,MAXBLOCKS)
  ! Plane average q^2 -> TKE = 1/2 * q^2  vs. z:
  real :: q2(NZB,MAXBLOCKS)
  ! Plane average of rms fluctuations vs. z:
  real :: urms2(NZB,MAXBLOCKS),vrms2(NZB,MAXBLOCKS),wrms2(NZB,MAXBLOCKS)
  ! Plane average of longitudinal and transverse derivatives vs. z:
  real :: dudx2(NZB,MAXBLOCKS),dudy2(NZB,MAXBLOCKS)
  ! Plane average of longitudinal and transverse microscales vs. z:
  real :: lambdaf(NZB,MAXBLOCKS),lambdag(NZB,MAXBLOCKS)
  ! Plane average of dissipation vs. z:
  real :: eps_nu(NZB,MAXBLOCKS)

  ! If LES
  ! Plane average of turbulent viscosity and SGS dissipation vs. z:
  real :: nut(NZB,MAXBLOCKS),eps_nut(NZB,MAXBLOCKS)

  ! Average presssure field, and rms fluctuation:
  real :: pavg(NZB,MAXBLOCKS),prms2(NZB,MAXBLOCKS)

  !,uavgaux,vavgaux,wavgaux,fconstn
  !real dissavg,dissavgaux

  real :: dux,dvy,dwz,duy,duz,dvx,dvz,dwx,dwy

  integer i,j,k,kk,lb,blockID,ierr

  logical, save :: initfile = .TRUE.

  character(80) :: filename



  !--- Initialize averages:
  zc(:,:)   = 0.; 
  uavg(:,:) = 0.; vavg(:,:) = 0.; wavg(:,:) = 0.;
  omgx2(:,:) = 0.; omgy2(:,:) = 0.; omgz2(:,:) = 0.;

  uu(:,:) = 0.; vv(:,:) = 0.; ww(:,:) = 0.; uv(:,:) = 0.; uw(:,:) = 0.; vw(:,:) = 0.;
  q2(:,:) = 0.; urms2(:,:) = 0.; vrms2(:,:) = 0.; wrms2(:,:) = 0.;

  dudx2(:,:) = 0.; dudy2(:,:) = 0.;
  lambdaf(:,:) = 0.; lambdag(:,:) = 0.;

  eps_nu(:,:) = 0.; 
  nut(:,:)    = 0.; eps_nut(:,:) = 0.;
  
  pavg(:,:) = 0.; prms2(:,:) = 0.;

  ! Call filling of cell centered velocities and vorticity:
  call ins_velomg2center(blockList,blockCount)


  ! Compute zc poisitions of cell centers:
  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz, coord, bsize:
     call Grid_getDeltas(blockID,del)
     call Grid_getBlkCenterCoords(blockID,coord)
     call Grid_getBlkPhysicalSize(blockID,bsize)
  
     dz = del(KAXIS)
     ! cell center z locations for the block:
     do k=1,NZB
        zc(k,lb) = coord(KAXIS) - 0.5*bsize(KAXIS) + real(k-1)*dz + 0.5*dz 
     enddo

  enddo

  ! For starters assume blockCount is 1, as in 
  ! uniform grid, nprocs in z direction: 
  ! -------------------------------------------------------------------
  ! Average velocities in cell centers of plane:
  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     ! Cell volume:
     volcell = del(IAXIS)*del(JAXIS)*del(KAXIS)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     ! Now, by plane go to cell centers and add velocities:
     kk =0
     do k=GRID_KLO,GRID_KHI

        kk = kk+1

        ! U average:
        uavg(kk,lb) = sum(solnData(VELX_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,k))*volcell

        ! V average:
        vavg(kk,lb) = sum(solnData(VELY_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,k))*volcell

        ! W average:
        wavg(kk,lb) = sum(solnData(VELZ_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,k))*volcell

        ! Pres average:
        pavg(kk,lb) = sum(solnData(PRES_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,k))*volcell

        do j=GRID_JLO,GRID_JHI
         do i=GRID_ILO,GRID_IHI

            ! Omgx2 average:
            omgx2(kk,lb) = omgx2(kk,lb) + solnData(OMGX_VAR,i,j,k)**2. * volcell                   

            ! Omgy2 average:
            omgy2(kk,lb) = omgy2(kk,lb) + solnData(OMGY_VAR,i,j,k)**2. * volcell

            ! Omgz2 average:
            omgz2(kk,lb) = omgz2(kk,lb) + solnData(OMGZ_VAR,i,j,k)**2. * volcell

            ! Velocity products:
            ! uu:
            uu(kk,lb) = uu(kk,lb) + solnData(VELX_VAR,i,j,k)*solnData(VELX_VAR,i,j,k)*volcell
            ! vv:
            vv(kk,lb) = vv(kk,lb) + solnData(VELY_VAR,i,j,k)*solnData(VELY_VAR,i,j,k)*volcell
            ! ww:
            ww(kk,lb) = ww(kk,lb) + solnData(VELZ_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k)*volcell
            ! uv:
            uv(kk,lb) = uv(kk,lb) + solnData(VELX_VAR,i,j,k)*solnData(VELY_VAR,i,j,k)*volcell
            ! uw:
            uw(kk,lb) = uw(kk,lb) + solnData(VELX_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k)*volcell 
            ! vw:
            vw(kk,lb) = vw(kk,lb) + solnData(VELY_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k)*volcell

            ! Pressure product:
            prms2(kk,lb) = prms2(kk,lb) + solnData(PRES_VAR,i,j,k)*solnData(PRES_VAR,i,j,k)*volcell

            ! Longitudinal derivatives:
            dux = (facexData(VELC_FACE_VAR,i+1,j,k) - facexData(VELC_FACE_VAR,i,j,k))/del(IAXIS)
            dvy = (faceyData(VELC_FACE_VAR,i,j+1,k) - faceyData(VELC_FACE_VAR,i,j,k))/del(JAXIS) 
            dwz = (facezData(VELC_FACE_VAR,i,j,k+1) - facezData(VELC_FACE_VAR,i,j,k))/del(KAXIS)

            dudx2(kk,lb) = dudx2(kk,lb) + (dux**2. + dvy**2. + dwz**2.)*volcell

            ! Transverse derivatives:
            duy = 0.25*(facexData(VELC_FACE_VAR,i+1,j+1,k) - facexData(VELC_FACE_VAR,i+1,j-1,k) + &
                        facexData(VELC_FACE_VAR,i,j+1,k)   - facexData(VELC_FACE_VAR,i,j-1,k) )/del(JAXIS)
            duz = 0.25*(facexData(VELC_FACE_VAR,i+1,j,k+1) - facexData(VELC_FACE_VAR,i+1,j,k-1) + &
                        facexData(VELC_FACE_VAR,i,j,k+1)   - facexData(VELC_FACE_VAR,i,j,k-1) )/del(KAXIS)

            dvx = 0.25*(faceyData(VELC_FACE_VAR,i+1,j+1,k) - faceyData(VELC_FACE_VAR,i-1,j+1,k) + &  
                        faceyData(VELC_FACE_VAR,i+1,j,k)   - faceyData(VELC_FACE_VAR,i-1,j,k) )/del(IAXIS)
            dvz = 0.25*(faceyData(VELC_FACE_VAR,i,j+1,k+1) - faceyData(VELC_FACE_VAR,i,j+1,k-1) + &
                        faceyData(VELC_FACE_VAR,i,j,k+1)   - faceyData(VELC_FACE_VAR,i,j,k-1) )/del(KAXIS)
            
            dwx = 0.25*(facezData(VELC_FACE_VAR,i+1,j,k+1) - facezData(VELC_FACE_VAR,i-1,j,k+1) + &
                        facezData(VELC_FACE_VAR,i+1,j,k)   - facezData(VELC_FACE_VAR,i-1,j,k) )/del(IAXIS)
            dwy = 0.25*(facezData(VELC_FACE_VAR,i,j+1,k+1) - facezData(VELC_FACE_VAR,i,j-1,k+1) + &
                        facezData(VELC_FACE_VAR,i,j+1,k)   - facezData(VELC_FACE_VAR,i,j-1,k) )/del(JAXIS)

            dudy2(kk,lb) = dudy2(kk,lb) + (duy**2. +duz**2. +dvx**2. +dvz**2. +dwx**2. +dwy**2.)*volcell 

         enddo
        enddo

        ! nut average:
        nut(kk,lb) = sum(solnData(TVIS_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,k))*volcell 

     enddo

     ! Sum to blocks Plane total volume:
     tvol(lb) = real(NXB*NYB*1)*volcell

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
           
  enddo
  
  ! Here information among processors and blocks containing the plane
  ! should be gathered

  ! Sum processors average velocities:
  !uavgaux = uavg
  !call MPI_Allreduce(uavgaux, uavg, 1, MPI_DOUBLE_PRECISION,&
  !                   MPI_SUM, MPI_COMM_WORLD, ierr) 
  !vavgaux = vavg
  !call MPI_Allreduce(vavgaux, vavg, 1, MPI_DOUBLE_PRECISION,&
  !                   MPI_SUM, MPI_COMM_WORLD, ierr) 
  !wavgaux = wavg
  !call MPI_Allreduce(wavgaux, wavg, 1, MPI_DOUBLE_PRECISION,&
  !                   MPI_SUM, MPI_COMM_WORLD, ierr) 
  ! Sum processors volumes
  !tvolaux = tvol  
  !call MPI_Allreduce(tvolaux, tvol, 1, MPI_DOUBLE_PRECISION,&
  !                   MPI_SUM, MPI_COMM_WORLD, ierr)


  ! Averages in the planes, per block:
  do lb = 1,blockCount

     ! Velocities:
     uavg(:,lb) = uavg(:,lb)/tvol(lb)
     vavg(:,lb) = vavg(:,lb)/tvol(lb)
     wavg(:,lb) = wavg(:,lb)/tvol(lb)

     ! Pressure:
     pavg(:,lb) = pavg(:,lb)/tvol(lb)

     ! Vorticity: 
     omgx2(:,lb) = omgx2(:,lb)/tvol(lb)
     omgy2(:,lb) = omgy2(:,lb)/tvol(lb)
     omgz2(:,lb) = omgz2(:,lb)/tvol(lb)

     ! Longitudinal and transverse derivatives dudx2 and dudy2
     dudx2(:,lb) = dudx2(:,lb)/tvol(lb)
     dudy2(:,lb) = dudy2(:,lb)/tvol(lb)

     ! Turbulent viscosity:
     nut(:,lb) = nut(:,lb)/tvol(lb)

  enddo

 
  ! Write up an average velocities file
  if (gr_meshMe .eq. gr_meshNumProcs/2) then

      if (initfile) then
         OPEN(UNIT=240,FILE='./IOData/AVGVELS_z.res',STATUS='REPLACE', &
              FORM='formatted')
         write(240,'(I8,4g14.6)')istep,time,sum(uavg(1:NZB,1))/real(NZB), &
                                            sum(vavg(1:NZB,1))/real(NZB), &
                                            sum(wavg(1:NZB,1))/real(NZB)
         close(240)
         initfile = .false. 
      else
         OPEN(UNIT=240,FILE='./IOData/AVGVELS_z.res',STATUS='OLD',     &
              FORM='formatted',POSITION='APPEND')
         write(240,'(I8,4g14.6)')istep,time,sum(uavg(1:NZB,1))/real(NZB), &
                                            sum(vavg(1:NZB,1))/real(NZB), &
                                            sum(wavg(1:NZB,1))/real(NZB)
         close(240)
      endif

   endif


  ! Keep computing:
  ! Now urms2 and q2:
  ! Substract Mean Velocities, pressure: 
  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)


     ! Reynolds Stresses:
     kk =0
     do k=GRID_KLO,GRID_KHI
        kk = kk+1

        uu(kk,lb) = uu(kk,lb)/tvol(lb) - uavg(kk,lb)*uavg(kk,lb)
        vv(kk,lb) = vv(kk,lb)/tvol(lb) - vavg(kk,lb)*vavg(kk,lb)
        ww(kk,lb) = ww(kk,lb)/tvol(lb) - wavg(kk,lb)*wavg(kk,lb)

        uv(kk,lb) = uv(kk,lb)/tvol(lb) - uavg(kk,lb)*vavg(kk,lb)
        uw(kk,lb) = uw(kk,lb)/tvol(lb) - uavg(kk,lb)*wavg(kk,lb)
        vw(kk,lb) = vw(kk,lb)/tvol(lb) - vavg(kk,lb)*wavg(kk,lb)

        prms2(kk,lb) = prms2(kk,lb)/tvol(lb) - pavg(kk,lb)*pavg(kk,lb)

     enddo

     ! Now, get velocity fluctuations:
     kk =0
     do k=GRID_KLO,GRID_KHI

        kk = kk+1
      
        ! u fluctuation:
        solnData(VELX_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,k) = &
        solnData(VELX_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,k) - uavg(kk,lb)       

        ! v fluctuation:
        solnData(VELY_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,k) = &
        solnData(VELY_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,k) - vavg(kk,lb)
        
        ! w fluctuation:
        solnData(VELZ_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,k) = &
        solnData(VELZ_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,k) - wavg(kk,lb)

        ! Rms2:
        do j=GRID_JLO,GRID_JHI
         do i=GRID_ILO,GRID_IHI

            ! urms2:
            urms2(kk,lb) = urms2(kk,lb) + solnData(VELX_VAR,i,j,k)*solnData(VELX_VAR,i,j,k)*volcell  

            ! vrms2:
            vrms2(kk,lb) = vrms2(kk,lb) + solnData(VELY_VAR,i,j,k)*solnData(VELY_VAR,i,j,k)*volcell 

            ! wrms2:
            wrms2(kk,lb) = wrms2(kk,lb) + solnData(VELZ_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k)*volcell

         enddo 
        enddo

        ! divide urms2, vrms2, wrms2 by tvol:
        urms2(kk,lb) = urms2(kk,lb)/tvol(lb)
        vrms2(kk,lb) = vrms2(kk,lb)/tvol(lb)
        wrms2(kk,lb) = wrms2(kk,lb)/tvol(lb)

        ! Get q^2:
        q2(kk,lb) = urms2(kk,lb) + vrms2(kk,lb) + wrms2(kk,lb)

     enddo

     ! Resolved dissipation:
     eps_nu(:,lb) = mvisc*(dudx2(:,lb) + dudy2(:,lb)) 


     ! Finally Taylor microscales:
     kk =0
     do k=GRID_KLO,GRID_KHI

        kk = kk+1

        ! Longitudinal Lf = sqrt(2/3 * q2 / (1/3 dudx2))
        lambdaf(kk,lb) = sqrt(2./3. * q2(kk,lb) / (1./3. * dudx2(kk,lb)))
     
        ! Transversal  Lg = sqrt(2/3 * q2 / (1/6 dudy2)) 
        lambdag(kk,lb) = sqrt(2./3. * q2(kk,lb) / (1./6. * dudy2(kk,lb)))

     enddo
   
     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
           
  enddo

  ! -------------------------------------------------------------------
  ! Write Processor file:
  write(filename,'("IOData/statsz.",i4.4,".",i4.4,".bin")') iwrite, gr_meshMe
  if (gr_meshMe .eq. MASTER_PE) then 
    write(*,*) ' '
    write(*,*) 'Writing filename=',trim(filename)
  endif
  open(unit=1,file=trim(filename),status='replace',form='unformatted')
  write(1) time
  write(1) istep,blockCount,NXB,NYB,NZB  
  do lb = 1,blockCount
     blockID = blockList(lb)
     
     ! Get blocks dx, dy ,dz, coord, bsize:
     call Grid_getDeltas(blockID,del)
     call Grid_getBlkCenterCoords(blockID,coord)
     call Grid_getBlkPhysicalSize(blockID,bsize)
     
     ! Write out deltas, coord, bsize:
     write(1) (del(i)  , i=1,MDIM)
     write(1) (coord(i), i=1,MDIM)
     write(1) (bsize(i), i=1,MDIM)

     ! Write slices zc location:
     write(1) (zc(kk,lb), kk=1,NZB)

     ! Write mean velocities, and pressure:
     write(1) (uavg(kk,lb), kk=1,NZB)
     write(1) (vavg(kk,lb), kk=1,NZB)
     write(1) (wavg(kk,lb), kk=1,NZB) 
     write(1) (pavg(kk,lb), kk=1,NZB)
 
     ! Write vorticities:    
     write(1) (omgx2(kk,lb), kk=1,NZB)
     write(1) (omgy2(kk,lb), kk=1,NZB)
     write(1) (omgz2(kk,lb), kk=1,NZB)

     ! Write q^2 and rms2 fluctuations:
     write(1) (q2(kk,lb), kk=1,NZB)
     write(1) (urms2(kk,lb), kk=1,NZB)
     write(1) (vrms2(kk,lb), kk=1,NZB)
     write(1) (wrms2(kk,lb), kk=1,NZB)
     write(1) (prms2(kk,lb), kk=1,NZB)

     ! Write dissipation and Taylor microscales:
     write(1) (eps_nu(kk,lb), kk=1,NZB)
     write(1) (lambdaf(kk,lb), kk=1,NZB)
     write(1) (lambdag(kk,lb), kk=1,NZB)

     ! Finally Reynolds stresses:
     write(1) (uu(kk,lb), kk=1,NZB)
     write(1) (vv(kk,lb), kk=1,NZB)
     write(1) (ww(kk,lb), kk=1,NZB)
     write(1) (uv(kk,lb), kk=1,NZB)
     write(1) (uw(kk,lb), kk=1,NZB)
     write(1) (vw(kk,lb), kk=1,NZB) 


  enddo

  close(1)

  return
  end subroutine


