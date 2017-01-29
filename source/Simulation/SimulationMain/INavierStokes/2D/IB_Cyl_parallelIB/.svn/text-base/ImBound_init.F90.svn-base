!!****if* source/Simulation/SimulationMain/INavierStokes/2D/IB_Cyl_parallelIB/ImBound_init
!!
!! NAME
!!
!!  ImBound_init
!!
!!
!! SYNOPSIS
!!
!!  call ImBound_init(restart)
!!  
!! VARIABLES
!!
!! restart = restart flag, logical.
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***

subroutine ImBound_init(restart)

  use Driver_data, ONLY: dr_simTime

  use Driver_interface, ONLY : Driver_abortFlash, &
       Driver_getMype, Driver_getComm

  use ImBound_data, only : ib_meshMe, ib_meshComm, ib_globalMe
  use ImBound_data
  use gr_sbData, only : gr_sbNumBodies, gr_sbBodyInfo

  use Grid_interface, only :  Grid_sbSelectMaster, &
    Grid_sbBroadcastParticles, Grid_getBoundboxCentroids

  use gr_sbInterface, ONLY: gr_sbCreateGroups, gr_sbCreateParticles

  implicit none
#include "constants.h"
#include "Flash.h"

  logical, INTENT(IN) :: restart

  ! Local variables:

  real :: xpt, ypt, dtheta, angle, tita, dsb
  integer :: i, b, ibd, nodelocpos, NumVertices, NumTriangles

  real :: L1,L2,n1(MDIM),n2(MDIM)


  call Driver_getMype(MESH_COMM,ib_meshMe)
  call Driver_getComm(MESH_COMM,ib_meshComm)
  call Driver_getMype(GLOBAL_COMM,ib_globalMe)

  tita  = omg*dr_simTime 

!#define DATAFLG 1

#ifdef DATAFLG 
!DATAFLG == OURS

  ib_nbd = 1

  ib_nela  = 628   !200 !380
  ib_nnoda = ib_nela + 1

  allocate(ib_ABODY(ib_nbd))
  allocate(xo(ib_nbd),yo(ib_nbd))  


  ibd = 1 
  ib_ABODY(ibd) % lb = 0
  ib_ABODY(ibd) % mb = ib_nnoda
  ib_ABODY(ibd) % elb = 0
  ib_ABODY(ibd) % emb = ib_nela
  ib_ABODY(ibd) % memflag = 1
  ib_ABODY(ibd) % ron = 0.

  xo(ibd) = 0. !0.5
  yo(ibd) = 0. !0.5

  ! Define positions and normals
  dtheta = 2.*PI/ib_nela

  angle = 0.
  ib_sb(1) = 0.
  do i=1,ib_nnoda


    ! Positions of markers, normals and arclength var
    ib_xbus(i) = Ro*cos(angle)
    ib_ybus(i) = Ro*sin(angle)

    ib_xb(i) = xo + ib_xbus(i)*cos(tita) - ib_ybus(i)*sin(tita)
    ib_yb(i) = yo + ib_xbus(i)*sin(tita) + ib_ybus(i)*cos(tita)

    ib_nxL(i) = ib_xb(i)/Ro
    ib_nyL(i) = ib_yb(i)/Ro

    if (i .gt. 1) then
       !dsb = sqrt( (ib_xb(i)-ib_xb(i-1))**2 + (ib_yb(i)-ib_yb(i-1))**2 )
       dsb = dtheta*Ro
       ib_sb(i) = ib_sb(i-1) + dsb

       ib_AELEM(1,i-1) = 2
       ib_AELEM(2,i-1) = i-1
       ib_AELEM(3,i-1) = i
    endif

    ! Set velocities and accelerations
    ib_ubd(i) = -omg*Ro*sin(angle) ! 0.
    ib_vbd(i) =  omg*Ro*cos(angle) ! 0.
    ib_ubdd(i)= 0.
    ib_vbdd(i)= 0.
 
    ! Update angle
    angle = angle + dtheta

  end do

!#define DATAFLG

#else 
!if defined DATAFLG 

  allocate(xo(gr_sbNumBodies),yo(gr_sbNumBodies))

 


  do b=1,gr_sbNumBodies

     write(*,*) 'b=',b,'NumBods=',gr_sbNumBodies

     NumVertices = gr_sbBodyInfo(b)%NumVertices

     xo(b) = 0.0
     yo(b) = 0.0

     gr_sbBodyInfo(b)%xo = xo(b)
     gr_sbBodyInfo(b)%yo = yo(b)


     !write(*,*) 'xo=',gr_sbBodyInfo(b)%xo,'yo=',gr_sbBodyInfo(b)%yo
     !pause

     gr_sbBodyInfo(b) % sb(1) = 0.
     do i=1,NumVertices

        gr_sbBodyInfo(b)%xb(i) = xo(b) + gr_sbBodyInfo(b)%xbus(i)*cos(tita) - gr_sbBodyInfo(b)%ybus(i)*sin(tita)
        gr_sbBodyInfo(b)%yb(i) = yo(b) + gr_sbBodyInfo(b)%xbus(i)*sin(tita) + gr_sbBodyInfo(b)%ybus(i)*cos(tita) + ao*sin(2.*PI*freq_t*dr_simTime) 

        ! Set velocities and accelerations
        gr_sbBodyInfo(b)%ubd(i) = -omg*gr_sbBodyInfo(b)%xbus(i)*sin(tita) - omg*gr_sbBodyInfo(b)%ybus(i)*cos(tita) ! 0.
        gr_sbBodyInfo(b)%vbd(i) =  omg*gr_sbBodyInfo(b)%xbus(i)*cos(tita) - omg*gr_sbBodyInfo(b)%ybus(i)*sin(tita) + (2.*PI*freq_t)*ao*cos(2.*PI*freq_t*dr_simTime) ! 0.
        gr_sbBodyInfo(b)%ubdd(i)= -omg**2.*gr_sbBodyInfo(b)%xbus(i)*cos(tita) + omg**2.*gr_sbBodyInfo(b)%ybus(i)*sin(tita)
        gr_sbBodyInfo(b)%vbdd(i)= -omg**2.*gr_sbBodyInfo(b)%xbus(i)*sin(tita) - omg**2.*gr_sbBodyInfo(b)%ybus(i)*cos(tita) - (2.*PI*freq_t)**2.*ao*sin(2.*PI*freq_t*dr_simTime)

       if (i .gt. 1) then
          dsb = sqrt( (gr_sbBodyInfo(b)%xb(i)-gr_sbBodyInfo(b)%xb(i-1))**2 + (gr_sbBodyInfo(b)%yb(i)-gr_sbBodyInfo(b)%yb(i-1))**2 )
          gr_sbBodyInfo(b)%sb(i) = gr_sbBodyInfo(b)%sb(i-1) + dsb
       endif


     enddo

     ! Normals:
     i = 1
     L1    = sqrt((gr_sbBodyInfo(b)%xb(NumVertices)-gr_sbBodyInfo(b)%xb(NumVertices-1))**2. + &
                  (gr_sbBodyInfo(b)%yb(NumVertices)-gr_sbBodyInfo(b)%yb(NumVertices-1))**2.)
     L2    = sqrt((gr_sbBodyInfo(b)%xb(i+1)-gr_sbBodyInfo(b)%xb(i))**2. + &
                  (gr_sbBodyInfo(b)%yb(i+1)-gr_sbBodyInfo(b)%yb(i))**2.) 

     n1(1) = L1**(-1.)*(gr_sbBodyInfo(b)%yb(NumVertices)-gr_sbBodyInfo(b)%yb(NumVertices-1))
     n1(2) =-L1**(-1.)*(gr_sbBodyInfo(b)%xb(NumVertices)-gr_sbBodyInfo(b)%xb(NumVertices-1))

     n2(1) = L2**(-1.)*(gr_sbBodyInfo(b)%yb(i+1)-gr_sbBodyInfo(b)%yb(i))
     n2(2) =-L2**(-1.)*(gr_sbBodyInfo(b)%xb(i+1)-gr_sbBodyInfo(b)%xb(i))

     gr_sbBodyInfo(b) % nxl(NumVertices) = 0.5*(n1(1)+n2(1))
     gr_sbBodyInfo(b) % nyl(NumVertices) = 0.5*(n1(2)+n2(2))

     gr_sbBodyInfo(b) % nxl(i) = 0.5*(n1(1)+n2(1))
     gr_sbBodyInfo(b) % nyl(i) = 0.5*(n1(2)+n2(2))

     do i =2,gr_sbBodyInfo(b)%NumVertices-1

       L1    = sqrt((gr_sbBodyInfo(b)%xb(i)-gr_sbBodyInfo(b)%xb(i-1))**2. + &
                    (gr_sbBodyInfo(b)%yb(i)-gr_sbBodyInfo(b)%yb(i-1))**2.)
       L2    = sqrt((gr_sbBodyInfo(b)%xb(i+1)-gr_sbBodyInfo(b)%xb(i))**2. + &
                    (gr_sbBodyInfo(b)%yb(i+1)-gr_sbBodyInfo(b)%yb(i))**2.) 

       n1(1) = L1**(-1.)*(gr_sbBodyInfo(b)%yb(i)-gr_sbBodyInfo(b)%yb(i-1))
       n1(2) =-L1**(-1.)*(gr_sbBodyInfo(b)%xb(i)-gr_sbBodyInfo(b)%xb(i-1))

       n2(1) = L2**(-1.)*(gr_sbBodyInfo(b)%yb(i+1)-gr_sbBodyInfo(b)%yb(i))
       n2(2) =-L2**(-1.)*(gr_sbBodyInfo(b)%xb(i+1)-gr_sbBodyInfo(b)%xb(i))

       gr_sbBodyInfo(b) % nxl(i) = 0.5*(n1(1)+n2(1))
       gr_sbBodyInfo(b) % nyl(i) = 0.5*(n1(2)+n2(2))
 
     enddo
     

  enddo


  call Grid_getBoundboxCentroids()
  call Grid_sbSelectMaster()
  call gr_sbCreateParticles()
  call Grid_sbBroadcastParticles()


  ! Write .dat file:
  if (ib_meshMe .eq. MASTER_PE) then
  open(unit=113,file='./IOData/CYL2.dat',form='formatted')
  write(113,'(I8)') gr_sbBodyInfo(1) % NumVertices
  write(113,'(I8)') gr_sbBodyInfo(1) % NumTriangles

  do i = 1,gr_sbBodyInfo(1) % NumVertices
      write(113,'(2E22.12)') gr_sbBodyInfo(1)%xbus(i),gr_sbBodyInfo(1)%ybus(i)
  enddo

  do i = 1,gr_sbBodyInfo(1) % NumTriangles
      write(113,'(2I8)') gr_sbBodyInfo(1)%AELEM(2,i),gr_sbBodyInfo(1)%AELEM(3,i)
  enddo
  close(113)
  endif

   !call Driver_abortFlash("Done with Sim init")

#endif





  ! Write out the geometry:
  if (ib_meshMe .eq. MASTER_PE) then

     open(unit=113,file='./IOData/geometry_cylinder.plt',form='formatted')
     write(113,'(A,G12.5,A)')'TEXT X=75,Y=5,F=HELV-BOLD,C=RED,T=" T = ',dr_simTime,'"'

     do ibd=1,ib_nbd

        nodelocpos = ib_ABODY(ibd)%lb


        !call int2char(ibd,index_ibd)
        write(113,'(A)') 'VARIABLES = "X" , "Y"'
        write(113,'(2(A,I8),A)')                 &
        'ZONE T=BodyCylinder, N=',              &
        ib_ABODY(ibd)%mb+1,', E=',ib_ABODY(ibd)%emb,        &
        ',DATAPACKING = POINT, ZONETYPE = FETRIANGLE'

        do i = ib_ABODY(ibd)%lb + 1,ib_ABODY(ibd)%lb + ib_ABODY(ibd)%mb
           write(113,'(2E16.8)') ib_xb(i),ib_yb(i)
        enddo
        xpt = xo(ibd); ypt = yo(ibd);
        write(113,'(2E16.8)') xpt,ypt
        write(113,'(A)') ' '
        Do i=ib_ABODY(ibd)%elb+1,ib_ABODY(ibd)%elb+ib_ABODY(ibd)%emb
           write(113,'(3(I8))')                        &
                  ib_AELEM(2,i)-nodelocpos,            &
                  ib_AELEM(3,i)-nodelocpos,ib_ABODY(ibd)%mb+1
        enddo

     enddo
     close(113)
    
  endif

  !call Driver_abortFlash("After writing geometry file..")

  !----

!!$  ! Allocate size of gr_sbBodyInterp
!!$  allocate(gr_sbBodyInterp(gr_sbNumBodies))


!!$  ! Allocate fields of gr_sbBodyInterp for MLS interpolation functions
!!$  LOOP_ALLOCATE : do b = 1,gr_sbNumBodies
!!$     
!!$     npart
!!$
!!$         real, allocatable, dimension(:) :: dsxu,dsyu
!!$         real, allocatable, dimension(:) :: dsxv,dsyv
!!$         real, allocatable, dimension(:) :: dsxu2,dsyv2 
!!$         integer,  allocatable, dimension(:,:) :: ielemu,jelemu
!!$         integer,  allocatable, dimension(:,:) :: ielemv,jelemv
!!$         real, allocatable, dimension(:) :: UL,VL
!!$         real, allocatable, dimension(:) :: FuL,FvL
!!$         real, allocatable, dimension(:) :: FuL2,indFuL,FvL2,indFvL
!!$         real, allocatable, dimension(:,:) :: phileu,philev
!!$         real, allocatable, dimension(:) :: xminb,xmaxb,yminb,ymaxb
!!$#if NDIM == 3
!!$         real, allocatable, dimension(:) :: dszu
!!$         real, allocatable, dimension(:) :: dszv
!!$         real, allocatable, dimension(:) :: dsxw,dsyw,dszw
!!$         integer, allocatable, dimension(:,:) :: kelemu
!!$         integer, allocatable, dimension(:,:) :: kelemv
!!$         integer, allocatable, dimension(:,:) :: kelemw
!!$         real, allocatable, dimension(:) :: WL
!!$         real, allocatable, dimension(:) :: FwL
!!$         real, allocatable, dimension(:) :: FwL2,indFwL
!!$         real, allocatable, dimension(:,:) :: philew
!!$         real, allocatable, dimension(:) :: zminb,zmaxb
!!$#endif
!!$
!!$
!!$  enddo LOOP_ALLOCATE


  if(ib_globalMe==MASTER_PE)print*,'Immersed_Boundaries initialized'
  return

end subroutine ImBound_init
