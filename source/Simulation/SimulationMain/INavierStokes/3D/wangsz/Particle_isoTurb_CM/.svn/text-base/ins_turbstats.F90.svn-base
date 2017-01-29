! subroutine ins_turbstats
! Computes turbulent statistics for isotropic turbulence simulations
!
! Marcos Vanella, July 2007
! --------------------------------------------------------------------------

  subroutine ins_turbstats(mype,istep,time,mvisc, &
                           blockCount,blockList,turbkin)

  ! Modules Use:
  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits

      
  implicit none
#include "constants.h"
#include "Flash.h"
#include "IncompNS.h"
  include 'mpif.h'

  integer :: mype, istep
  real :: time, mvisc, turbkin
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(blockCount) :: blockList

  integer, parameter ::  ng = NGUARD
  integer, parameter ::  nxi= NGUARD + NXB
  integer, parameter ::  nyj= NGUARD + NYB
  integer, parameter ::  nzk= NGUARD + NZB
  integer, parameter ::  nxc= NGUARD + NXB + 1
  integer, parameter ::  nyc= NGUARD + NYB + 1
  integer, parameter ::  nzc= NGUARD + NZB + 1

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
  
  real del(3),dx,dy,dz
  real turbkinaux,volcell,tvol,tvolaux
  real varl,varlaux,vart,vartaux,dudx_v,dudy_v
  real lambdaf,lambdag
  real uavg,vavg,wavg,uavgaux,vavgaux,wavgaux,fconstn
  real dissavg,dissavgaux

  integer lb,blockID,ierr

  real, parameter :: diss = 0.4

  logical, save :: initfile = .TRUE.

  ! -------------------------------------------------------------------
  ! Average velocities, Check there is no drift:
  uavg = 0.
  vavg = 0.
  wavg = 0.
  tvol = 0.

  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     ! Cell volume:
     volcell = del(DIR_X)*del(DIR_Y)*del(DIR_Z) ! dx*dy*dz

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     ! U average:
     uavg = uavg + &
      sum(facexData(VELC_FACE_VAR,ng+2:nxc,ng+1:nyj,ng+1:nzk))*volcell

     ! V average:
     vavg = vavg + &
      sum(faceyData(VELC_FACE_VAR,ng+1:nxi,ng+2:nyc,ng+1:nzk))*volcell

     ! W average:
     wavg = wavg + &
      sum(facezData(VELC_FACE_VAR,ng+1:nxi,ng+1:nyj,ng+2:nzc))*volcell

     ! Sum to processors total volume:
     tvol = tvol + real(NXB*NYB*NZB)*volcell

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
           
  enddo
  

  ! Sum processors average velocities:
  uavgaux = uavg
  call MPI_Allreduce(uavgaux, uavg, 1, MPI_DOUBLE_PRECISION,&
                     MPI_SUM, MPI_COMM_WORLD, ierr) 

  vavgaux = vavg
  call MPI_Allreduce(vavgaux, vavg, 1, MPI_DOUBLE_PRECISION,&
                     MPI_SUM, MPI_COMM_WORLD, ierr) 

  wavgaux = wavg
  call MPI_Allreduce(wavgaux, wavg, 1, MPI_DOUBLE_PRECISION,&
                     MPI_SUM, MPI_COMM_WORLD, ierr) 


  ! Sum processors volumes
  tvolaux = tvol  
  call MPI_Allreduce(tvolaux, tvol, 1, MPI_DOUBLE_PRECISION,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)


  ! Average velocities:
  uavg = uavg/tvol
  vavg = vavg/tvol
  wavg = wavg/tvol


  if (mype .eq. 0) then

      if (initfile) then
         OPEN(UNIT=240,FILE='./IOData/AVGVELS.res',STATUS='REPLACE', &
              FORM='formatted')
         write(240,'(I8,4g14.6)')istep,time,uavg,vavg,wavg
         close(240)
      
      else
         OPEN(UNIT=240,FILE='./IOData/AVGVELS.res',STATUS='OLD',     &
              FORM='formatted',POSITION='APPEND')
         write(240,'(I8,4g14.6)')istep,time,uavg,vavg,wavg
         close(240)
      endif

   endif

  ! Substract Mean Velocities: 
  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
     
     facexData(VELC_FACE_VAR,:,:,:) = facexData(VELC_FACE_VAR,:,:,:) - uavg
     faceyData(VELC_FACE_VAR,:,:,:) = faceyData(VELC_FACE_VAR,:,:,:) - vavg
     facezData(VELC_FACE_VAR,:,:,:) = facezData(VELC_FACE_VAR,:,:,:) - wavg

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
           
  enddo


  ! -------------------------------------------------------------------


  ! -------------------------------------------------------------------
  ! Compute Turbulent Kinetic Energy and Taylor Microscales:
  turbkin = 0.
  tvol    = 0.
  varl    = 0.
  vart    = 0.

  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     dx = del(DIR_X)
     dy = del(DIR_Y)
     dz = del(DIR_Z)

     ! Cell volume:
     volcell = dx*dy*dz

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)


     ! Turbulent kinnetic energy in cell centers of block,
     ! zero mean velocities:
!!$     unk(3,ng+1:nxi,ng+1:nyj,ng+1:nzk,lb) = .5/4.*(         &
!!$          (facevarx(1,ng+1:nxi,ng+1:nyj,ng+1:nzk,lb)     +  &
!!$           facevarx(1,ng+2:nxc,ng+1:nyj,ng+1:nzk,lb))**2 +  &
!!$          (facevary(1,ng+1:nxi,ng+1:nyj,ng+1:nzk,lb)     +  &
!!$           facevary(1,ng+1:nxi,ng+2:nyc,ng+1:nzk,lb))**2 +  &
!!$          (facevarz(1,ng+1:nxi,ng+1:nyj,ng+1:nzk,lb)     +  &
!!$           facevarz(1,ng+1:nxi,ng+1:nyj,ng+2:nzc,lb))**2 )               

      solnData(DUST_VAR,ng+1:nxi,ng+1:nyj,ng+1:nzk) = .5*(            & 
           (facexData(VELC_FACE_VAR,ng+2:nxc,ng+1:nyj,ng+1:nzk))**2 + &
           (faceyData(VELC_FACE_VAR,ng+1:nxi,ng+2:nyc,ng+1:nzk))**2 + &
           (facezData(VELC_FACE_VAR,ng+1:nxi,ng+1:nyj,ng+2:nzc))**2 ) 

      ! Sum to turbkin, cell value * cell volume:
      turbkin = turbkin + & 
         sum(solnData(DUST_VAR,ng+1:nxi,ng+1:nyj,ng+1:nzk))*volcell


      ! Longitudinal derivatives squared * cell volume: 
      solnData(DUST_VAR,ng+1:nxi,ng+1:nyj,ng+1:nzk) = 1.*(                        &
              (dx**(-1.)*(facexData(VELC_FACE_VAR,ng+2:nxc,ng+1:nyj,ng+1:nzk) -     &
                        facexData(VELC_FACE_VAR,ng+1:nxi,ng+1:nyj,ng+1:nzk)))**2  &
           +  (dy**(-1.)*(faceyData(VELC_FACE_VAR,ng+1:nxi,ng+2:nyc,ng+1:nzk) -     &
                        faceyData(VELC_FACE_VAR,ng+1:nxi,ng+1:nyj,ng+1:nzk)))**2  &
           +  (dz**(-1.)*(facezData(VELC_FACE_VAR,ng+1:nxi,ng+1:nyj,ng+2:nzc) -     &
                        facezData(VELC_FACE_VAR,ng+1:nxi,ng+1:nyj,ng+1:nzk)))**2)


      varl  = varl + &
              sum(solnData(DUST_VAR,ng+1:nxi,ng+1:nyj,ng+1:nzk))*volcell


      ! Transverse derivatives squared * cell volume:
      solnData(DUST_VAR,ng+1:nxi,ng+1:nyj,ng+1:nzk) = .5*(                        &
              (dy**(-1.)*(facexData(VELC_FACE_VAR,ng+2:nxc,ng+2:nyc,ng+1:nzk) -      &  ! du/dy**2
                       facexData(VELC_FACE_VAR,ng+2:nxc,ng+1:nyj,ng+1:nzk)))**2   &
           +  (dz**(-1.)*(facexData(VELC_FACE_VAR,ng+2:nxc,ng+1:nyj,ng+2:nzc) -      &  ! du/dz**2
                       facexData(VELC_FACE_VAR,ng+2:nxc,ng+1:nyj,ng+1:nzk)))**2   &           
           +  (dx**(-1.)*(faceyData(VELC_FACE_VAR,ng+2:nxc,ng+2:nyc,ng+1:nzk) -      &  ! dv/dx**2
                       faceyData(VELC_FACE_VAR,ng+1:nxi,ng+2:nyc,ng+1:nzk)))**2   &
           +  (dz**(-1.)*(faceyData(VELC_FACE_VAR,ng+1:nxi,ng+2:nyc,ng+2:nzc) -      &  ! dv/dz**2
                       faceyData(VELC_FACE_VAR,ng+1:nxi,ng+2:nyc,ng+1:nzk)))**2   &
           +  (dx**(-1.)*(facezData(VELC_FACE_VAR,ng+2:nxc,ng+1:nyj,ng+2:nzc) -      &  ! dw/dx**2
                       facezData(VELC_FACE_VAR,ng+1:nxi,ng+1:nyj,ng+2:nzc)))**2   &
           +  (dy**(-1.)*(facezData(VELC_FACE_VAR,ng+1:nxi,ng+2:nyc,ng+2:nzc) -      &  ! dw/dy**2
                       facezData(VELC_FACE_VAR,ng+1:nxi,ng+1:nyj,ng+2:nzc)))**2)                

      vart  = vart + &
              sum(solnData(DUST_VAR,ng+1:nxi,ng+1:nyj,ng+1:nzk))*volcell

      ! Sum to processors total volume:
      tvol = tvol + real(NXB*NYB*NZB)*volcell

      ! Release pointers:
      call Grid_releaseBlkPtr(blockID,solnData,CENTER)
      call Grid_releaseBlkPtr(blockID,facexData,FACEX)
      call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
      call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  enddo

  ! Sum to processors total dissipation dissavg = nu * <(dui/dxj)^2>:
  dissavg = mvisc*(varl + 2. * vart)

  ! Sum processors turbulent dissipation
  dissavgaux = dissavg
  call MPI_Allreduce(dissavgaux, dissavg, 1, MPI_DOUBLE_PRECISION,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)  
       
  ! Sum processors turbulent kinetic energies (q^2):
  turbkinaux = turbkin
  call MPI_Allreduce(turbkinaux, turbkin, 1, MPI_DOUBLE_PRECISION,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)              

  ! Sum processors longitudinal derivatives squared sum:
  varlaux = varl
  call MPI_Allreduce(varlaux, varl, 1, MPI_DOUBLE_PRECISION,      &
                       MPI_SUM, MPI_COMM_WORLD, ierr)              

  ! Sum processors transversal derivatives squared sum:
  vartaux = vart
  call MPI_Allreduce(vartaux, vart, 1, MPI_DOUBLE_PRECISION,      &
                     MPI_SUM, MPI_COMM_WORLD, ierr)    

  ! Sum processors volumes
  tvolaux = tvol  
  call MPI_Allreduce(tvolaux, tvol, 1, MPI_DOUBLE_PRECISION,      &
                     MPI_SUM, MPI_COMM_WORLD, ierr)

  ! Volume Averaged turbulent kinetic energy:
  turbkin = turbkin/tvol

  ! Volume Averaged longitudinal derivatives squared:
  dudx_v = varl/(3.*tvol)

  ! Volume Averaged transverse derivatives squared:
  dudy_v = vart/(3.*tvol)

  ! Longitudinal Taylor Microscale:
  lambdaf = sqrt(4./3. * turbkin / dudx_v)  ! Uses q^2 = 2*turbkin

  ! Transverse Taylor Microscale:
  lambdag = sqrt(4./3. * turbkin / dudy_v)  ! Uses q^2 = 2*turbkin

  ! Dissipation eps
  dissavg = dissavg/tvol 

  !!!! Remember, Eddy turnover time = turbkin / dissavg !!!!
  if (mype .eq. 0) then

     if (initfile) then

      OPEN(UNIT=241,FILE='./IOData/TURBKINLAM.res',STATUS='REPLACE', &
           FORM='formatted')
      write(241,'(I8,5g14.6)')istep,time,turbkin,lambdaf,lambdag,dissavg
      write(*,*) &
      'ISTEP      TIME      TURBKIN        Lambdaf       Lambdag        dissavg'
      write(*,'(I6,5g14.6)')istep,time,turbkin,lambdaf,lambdag,dissavg
      close(241)      

      initfile = .FALSE.
      
      else

      OPEN(UNIT=241,FILE='./IOData/TURBKINLAM.res',STATUS='OLD',&
           FORM='formatted',POSITION='APPEND')
      write(241,'(I8,5g14.6)')istep,time,turbkin,lambdaf,lambdag,dissavg
      write(*,*) &
      'ISTEP      TIME      TURBKIN        Lambdaf       Lambdag        dissavg'
      write(*,'(I6,5g14.6)')istep,time,turbkin,lambdaf,lambdag,dissavg
      close(241) 

     endif

  endif
  ! -------------------------------------------------------------------

  ! Add Mean Velocities: 
  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
     
     facexData(VELC_FACE_VAR,:,:,:) = facexData(VELC_FACE_VAR,:,:,:) + uavg
     faceyData(VELC_FACE_VAR,:,:,:) = faceyData(VELC_FACE_VAR,:,:,:) + vavg
     facezData(VELC_FACE_VAR,:,:,:) = facezData(VELC_FACE_VAR,:,:,:) + wavg

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
           
  enddo

  return
  end subroutine


