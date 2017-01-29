!!****if* source/Simulation/SimulationMain/INavierStokes/2D/IB_Cyl/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! ARGUMENTS
!!
!!   
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for INS-isotropic turbulence problem.
!!
!!***

subroutine Simulation_init()

  use Grid_data, only : gr_meshMe

  use Driver_data, ONLY: dr_simTime

  use Driver_interface, ONLY : Driver_abortFlash

  use Simulation_data, ONLY : sim_xMin, sim_yMin, &
                              sim_xMax, sim_yMax, sim_gCell

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use ImBound_data

  implicit none

#include "constants.h"
#include "Flash.h"

  
  real :: xpt, ypt, dtheta, angle, tita, dsb
  integer :: i, ibd, nodelocpos

!  real :: omg, freq_nat, freq_t, ao

  call RuntimeParameters_get('xmin',    sim_xMin)
  call RuntimeParameters_get('ymin',    sim_yMin)
  call RuntimeParameters_get('xmax',    sim_xMax)
  call RuntimeParameters_get('ymax',    sim_yMax)

  sim_gCell = .true.


  ! Cylinder setup IB variables:
  Ro = 0.5 !0.5
  omg= 0.0 !2.0
  freq_nat = 0.*0.196 ! Natural shedding frequency at Re=100
  freq_t   = 1.0*freq_nat ! Translational motion frequency.
  ao = 0.2*2.*Ro ! Translational motion amplitude.


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

  xo(ibd) = 0.0 !0.5
  yo(ibd) = 0.0 !0.5

  ! Define positions and normals
  dtheta = 2.*PI/ib_nela

  angle = 0.
  tita  = omg*dr_simTime 
  ib_sb(1) = 0.
  do ibd=1,ib_nbd
  do i=ib_ABODY(ibd)%lb+1,ib_ABODY(ibd)%lb+ib_ABODY(ibd)%mb


    ! Positions of markers, normals and arclength var
    ib_xbus(i) = Ro*cos(angle)
    ib_ybus(i) = Ro*sin(angle)

    ib_xb(i) = xo(ibd) + ib_xbus(i)*cos(tita) - ib_ybus(i)*sin(tita)
    ib_yb(i) = yo(ibd) + ib_xbus(i)*sin(tita) + ib_ybus(i)*cos(tita) + ao*sin(2.*PI*freq_t*dr_simTime)

    ib_nxL(i) = ib_xb(i)/Ro
    ib_nyL(i) = ib_yb(i)/Ro

    if (i .gt. 1) then
       dsb = sqrt( (ib_xb(i)-ib_xb(i-1))**2 + (ib_yb(i)-ib_yb(i-1))**2 )
       !dsb = dtheta*Ro
       ib_sb(i) = ib_sb(i-1) + dsb

       ib_AELEM(1,i-1) = 2
       ib_AELEM(2,i-1) = i-1
       ib_AELEM(3,i-1) = i
    endif

    ! Set velocities and accelerations
    ib_ubd(i) = -omg*ib_xbus(i)*sin(tita) - omg*ib_ybus(i)*cos(tita) ! 0.
    ib_vbd(i) =  omg*ib_xbus(i)*cos(tita) - omg*ib_ybus(i)*sin(tita) + (2.*PI*freq_t)*ao*cos(2.*PI*freq_t*dr_simTime) ! 0.
    ib_ubdd(i)= -omg**2.*ib_xbus(i)*cos(tita) + omg**2.*ib_ybus(i)*sin(tita)
    ib_vbdd(i)= -omg**2.*ib_xbus(i)*sin(tita) - omg**2.*ib_ybus(i)*cos(tita) - (2.*PI*freq_t)**2.*ao*sin(2.*PI*freq_t*dr_simTime)
 
    ! Update angle
    angle = angle + dtheta

  end do
  end do

  ! Write .dat file:
  if (gr_meshMe .eq. MASTER_PE) then
  open(unit=113,file='./IOData/CYL.dat',form='formatted')
  write(113,'(I8)') ib_ABODY(1) % mb
  write(113,'(I8)') ib_ABODY(1) % emb

  do i = ib_ABODY(1)%lb+1,ib_ABODY(1)%lb+ib_ABODY(1)%mb
      write(113,'(2E22.12)') ib_xbus(i),ib_ybus(i)
  enddo

  do i = ib_ABODY(1)%elb+1,ib_ABODY(1)%elb+ib_ABODY(1)%emb
      write(113,'(2I8)') ib_AELEM(2,i),ib_AELEM(3,i)
  enddo
  close(113)
  endif

!  call Driver_abortFlash("Done with Sim init")

  ! Write out the geometry:
  if (gr_meshMe .eq. MASTER_PE) then

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



end subroutine Simulation_init
