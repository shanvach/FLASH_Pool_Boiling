#include "Flash.h"
#include "constants.h"

subroutine sim_writeOutputGrid
  use sim_outputGridData
  implicit none

  integer :: i,j,k, v, ov, outLUnit, nvarsOgWrite

  real,allocatable :: xOgCenter(:),yOgCenter(:),zOgCenter(:)
  real,allocatable :: xOgFaces(:),yOgFaces(:),zOgFaces(:)

  character(len=4) :: unkNames(NUNK_VARS)

  integer, parameter, dimension(1:nvarsOgOut) :: outvars = (/&
       NUP1_VAR,TEMP_VAR,PRES_VAR,DENS_VAR,VELX_VAR,VELZ_VAR,SUMX_MSCALAR,YE_MSCALAR,SUMY_MSCALAR,&
       NEU_SPEC,H_SPEC,HE_SPEC,LI_SPEC,BE_SPEC, B_SPEC, C_SPEC, N_SPEC, O_SPEC,&
       F_SPEC,NE_SPEC,NA_SPEC,MG_SPEC,AL_SPEC,SI_SPEC, P_SPEC, S_SPEC,&
       CL_SPEC,AR_SPEC, K_SPEC,CA_SPEC,SC_SPEC,TI_SPEC, V_SPEC,CR_SPEC,&
       MN_SPEC,FE_SPEC,CO_SPEC,NI_SPEC,&
       CU_SPEC,ZN_SPEC,GA_SPEC,GE_SPEC,&
       NI56_MSCALAR, NUMP_VAR, NUP0_VAR,NUMC_VAR/)


  do v=1,NUNK_VARS
     call Simulation_mapIntToStr(v, unkNames(v), MAPBLOCK_UNK)
     
  end do


  do k=kOgB0,kOgE0
     do j=jOgB0,jOgE0
        do i=0,nIOg+1
            outGridData(i,j,k,NUMC_VAR) = outGridNumCount(i,j,k)

            if (outGridNumCount(i,j,k) .NE. 0) then
               do v=1,nvarsOg
                  if ((v .NE. NUMC_VAR) .AND. (v .NE. NUP1_VAR)) then
                     outGridData(i,j,k,v) = outGridData(i,j,k,v) / outGridMappedVol(i,j,k)
                  end if
               end do
            end if
        end do
     end do
  end do


  allocate(xOgCenter(0:nIOg+1))
  allocate(yOgCenter(0:nJOg+1))
  allocate(zOgCenter(0:nKOg+1))
  allocate(xOgFaces(0:nIOg+2))
  allocate(yOgFaces(0:nJOg+2))
  allocate(zOgFaces(0:nKOg+2))

  if (xOgMin .LE. 0.0) then     !linear spacing...

     xOgStep = (xOgMax - xOgMin) / nIOg
     do i=0,nIOg+1
        xOgCenter(i) = xOgMin + (real(i)-0.5) * xOgStep 
     end do
     do i=0,nIOg+2
        xOgFaces(i) = xOgMin + (real(i)-1.0) * xOgStep 
     end do

  else                          !log spacing...
     xOgFact = 10.0**((alog10(xOgMax) - alog10(xOgMin)) / nIOg)
     do i=0,nIOg+1
        xOgCenter(i) = xOgMin * xOgFact**(real(i)-0.5)
     end do
     do i=0,nIOg+2
        xOgFaces(i) = xOgMin * xOgFact**(real(i)-1.0)
     end do
     do i=0,nIOg+1
        xOgCenter(i) = 0.5*(xOgFaces(i)+xOgFaces(i+1))
     end do
  end if

  print*,'ndimOg=',ndimOg
  if (ndimOg>1) then
     yOgStep = (yOgMax - yOgMin) / nJOg
     print*,'yOgMin=',yOgMin
     print*,'yOgMax=',yOgMax
     print*,'yOgStep=',yOgStep
     print*,'jOgB0,jOgE0=',jOgB0,jOgE0

     do i=jOgB0,jOgE0
        yOgFaces(i) = yOgMin + (i-1) * yOgStep 
        yOgCenter(i) = yOgMin + (i-0.5) * yOgStep 
        print*,i,yOgFaces(i),yOgCenter(i)
     end do
     yOgFaces(nJOg+1) = yOgMax
     yOgFaces(jOgE0+1) = yOgMax + yOgStep
  end if

  if (ndimOg>2) then
     zOgStep = (zOgMax - zOgMin) / nKog
     do i=kOgB0,kOgE0
        zOgFaces(i) = zOgMin + (i-1) * zOgStep 
        zOgCenter(i) = zOgMin + (i-0.5) * zOgStep 
     end do
     zOgFaces(nKog+1) = zOgMax
     zOgFaces(kOgE0+1) = zOgMax + zOgStep
  end if


990 format("FLASH-to-radTran FMT v003.002 2009-10-28a 1D")
20990 format("FLASH-to-radTran FMT v003.002 2009-10-28a 2D")
30990 format("FLASH-to-radTran FMT v003.001 2009-10-28a 3D")
99430 format("#","  I ","  J ","  K ",100(2x,a4,3x))
99130 format("#","    ","    ","    ",100(I6,3x))
!!994 format("#","  I ","       R       ",9(3x,a5,3x),100(1x,a5,1x))
!!991 format("#","    ","               ",9(I8,3x),100(I4,3x))
!!!994 format("#","  I ","       R       ",5(3x,a5,3x),100(1x,a5,1x))
!!!991 format("#","    ","               ",5(I8,3x),100(I4,3x))
994 format("#","  I ","     R     ","    R_left     ","    R_right    ",5(3x,a5,3x),2(4x,a6,4x),100(4x,a5,2x))
991 format("#","    ","           ","               ","               ",5(I8,3x),2(I11,3x),100(I8,3x))
20994 format("#","  I ","  J ","     R     ","     Z     ",&
           "    R_left     ","    R_right    ","    Z_left     ","    Z_right    ",6(3x,a5,3x),2(4x,a6,4x),100(4x,a5,2x))
30994 format("#","  I ","  J ","  K ","     X     ","     Y     ","     Z     ",&
           "    X_left     ","    X_right    ","    Y_left     ","    Y_right    ","    Z_left     ","    Z_right    ",&
           5(3x,a5,3x),2(4x,a6,4x),100(4x,a5,2x))
99230 format(1x,3(I4),100(G9.2))
!!992 format(1x,I4,G15.5,9(1x,1PG10.3),0P,100(1x,F6.2))
!!!992 format(1x,I4,1PG15.5,5(1x,1PG10.3),0P,100(1x,F6.2))
992 format(1x,I4,1x,1PG10.3,1PG15.5,1PG15.5,5(1x,1PG10.3),2(1x,0PF13.10),1P,100(1x,G10.3))
20992 format(1x,2I4,2(1x,1PG10.3),2(1PG15.5,1PG15.5),6(1x,1PG10.3),2(1x,0PF13.10),1P,100(1x,G10.3))
30992 format(1x,3I4,3(1x,1PG10.3),3(1PG15.5,1PG15.5),5(1x,1PG10.3),2(1x,0PF13.10),1P,100(1x,G10.3))

  nvarsOgWrite = nvarsOgOut ! min(nvarsOg,19)
  outLUnit = 80
  open(outLUnit,file=sim_ogRadTranDataFileName,status='UNKNOWN')
  select case (ndimOg)
  case(1)
     write(outLUnit,990)
     write(outLUnit,994) (unkNames(outVars(v)),v=1,nvarsOgWrite)
  case(2)
     write(outLUnit,20990)
     write(outLUnit,20994) (unkNames(outVars(v)),v=1,nvarsOgWrite)
  case(3)
     write(outLUnit,30990)
     write(outLUnit,30994) (unkNames(outVars(v)),v=1,nvarsOgWrite)
  end select
!!$  write(outLUnit,991) (v,v=1,nvarsOgWrite)
  do k=kOgB0,kOgE0
     do j=jOgB0,jOgE0
        do i=0,nIOg+1
!           if (outGridData(i,j,k,NUMP_VAR) .NE. 0.0) then
              do v=1,nvarsOgWrite
                 ov = outVars(v)
                    if ((ov == NI56_MSCALAR .OR. ov == SUMY_MSCALAR .OR. &
                      (ov .GE. SPECIES_BEGIN .AND. ov .LE. SPECIES_END))) then
                    if (outGridData(i,j,k,ov) .GE. 0 .AND. &
                         outGridData(i,j,k,ov) .LT. 1.0E-99) then
                       outLineData(v) = 0.0
                    else if (outGridData(i,j,k,ov) .LT. 0) then

8888                   format('WARNING negative data in "',A,'" at',I3,',',I3,',',I3,':',G30.20)
                       print 8888,unkNames(ov),i,j,k,outGridData(i,j,k,ov)
                       outLineData(v) = 0.0
                    else
                       outLineData(v) = outGridData(i,j,k,ov)
                    end if
                 else
                    outLineData(v) = outGridData(i,j,k,ov)
                 end if
              end do
              select case (ndimOg)
              case(1)
                 write(outLUnit,992) i,xOgCenter(i),xOgFaces(i),xOgFaces(i+1),(outLineData(v),v=1,nvarsOgWrite)
              case(2)
                 write(outLUnit,20992) i,j,xOgCenter(i),yOgCenter(j), &
                      xOgFaces(i),xOgFaces(i+1), &
                      yOgFaces(j),yOgFaces(j+1), &
                      (outLineData(v),v=1,nvarsOgWrite)
              case(3)
                 write(outLUnit,30992) i,j,k,xOgCenter(i),yOgCenter(j),zOgCenter(k), &
                      xOgFaces(i),xOgFaces(i+1), &
                      yOgFaces(j),yOgFaces(j+1), &
                      zOgFaces(k),zOgFaces(k+1), &
                      (outLineData(v),v=1,nvarsOgWrite)
              end select
!           end if
        end do
     end do
  end do

  deallocate(xOgCenter)
  deallocate(yOgCenter)
  deallocate(zOgCenter)
  deallocate(xOgFaces)
  deallocate(yOgFaces)
  deallocate(zOgFaces)



end subroutine sim_writeOutputGrid
