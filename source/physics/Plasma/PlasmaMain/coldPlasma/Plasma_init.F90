subroutine Plasma_init(blockCount,blockList,restart)

! Read runtime parameters and initialize constants declared in Plasma_data
! script

   use Plasma_data
   use RuntimeParameters_interface, ONLY : RuntimeParameters_get
   use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs, &
                                Driver_getComm, Driver_getNstep

   use Grid_interface, only: Grid_getDeltas, Grid_getBlkIndexLimits,&
                             Grid_getBlkPtr, Grid_releaseBlkPtr,    &
                             Grid_fillGuardCells

   use Plasma_interface, only: Plasma_Feed, Plasma_sumIons,&
                               Plasma_getNorm, Plasma_netCharge, Plasma_elPotential

   implicit none

  include 'Flash_mpi.h'
#include "constants.h"
#include "Flash.h"
#include "Plasma.h"

   integer, INTENT(INOUT) :: blockCount
   integer, INTENT(INOUT) :: blockList(MAXBLOCKS)
   logical, INTENT(IN)    :: restart

   ! Local Variables
   integer ::  blockID,lb,i,j
   real ::  del(MDIM)
   integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
   logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
   real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

   real, dimension(1) :: rand_noise
   real :: plasma_source !scalar, plasma source rate for species m-3 s-1
   real :: nrel, nrna, nrni !nrh0, nrh1, nrh2, nrh6, nrh7, nrh8, nrh9 

   real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: oldPhi
   real :: tempDT, Troom, rrc, rrx, rry, trc, vrx, vry, pu(20), pt(20)

   integer :: fac, pcount

   fac    = 9
   tempDT = 1.0
   Troom  = 0.0

   pu(1) = -0.0000
   pu(2) =  0.0005
   pu(3) = -0.0035
   pu(4) =  0.0065
   pu(5) =  0.0267
   pu(6) = -0.1211
   pu(7) =  0.0232
   pu(8) =  0.4247
   pu(9) =  0.0003

   pt(1) =  0.0000
   pt(2) = -0.0002
   pt(3) =  0.0021
   pt(4) = -0.0106
   pt(5) =  0.0281
   pt(6) = -0.0245
   pt(7) = -0.0435
   pt(8) = -0.0030
   pt(9) =  0.5303

   !pu(1) = -0.0000
   !pu(2) =  0.0011
   !pu(3) = -0.0135
   !pu(4) =  0.0743
   !pu(5) = -0.1839
   !pu(6) =  0.0650
   !pu(7) =  0.4133
   !pu(8) =  0.0009
  
   !pt(1) = -0.0000
   !pt(2) =  0.0004
   !pt(3) = -0.0032
   !pt(4) =  0.0102
   !pt(5) = -0.0010
   !pt(6) = -0.0592
   !pt(7) =  0.0012
   !pt(8) =  0.5300
 
   nrel = 0.99*1e18       ! Electrons
   nrna = 1e26            ! neutral
   nrni = 1e18            ! ions

   call Driver_getMype(MESH_COMM, pls_meshMe)
   call Driver_getNumProcs(MESH_COMM, pls_meshNumProcs)
   call Driver_getComm(MESH_COMM, pls_meshComm)

   call RuntimeParameters_get("cflflg", pls_cflflg)
   call RuntimeParameters_get("cfl", pls_cfl)
   call RuntimeParameters_get("sigma",pls_sigma)
   call RuntimeParameters_get("dtspec",pls_dtspec)
   call RuntimeParameters_get("vel_prolong_method",pls_prol_method)
   call RuntimeParameters_get("pct_noise",pls_pct_noise)

   call Driver_getNstep(pls_nstep)
   pls_restart=restart

   pls_dcoeff = 1e3 ! test electron diffusion coefficient

   ! He+, N2+, O-, O2+
   pls_Cmi_net = 6.6464764e-27 + 2*2.3258671e-26 + &
                 2.6566962e-26 + 2*2.6566962e-26 
   
   pls_epsilon0 = 8.85418782e-12
   pls_l2target = 1e-7

   pls_Ckb = 1.38064852e-23
   pls_Cme = 9.10938356e-31
   pls_Ce  = 1.60217662e-19
      pls_gam = EXP(0.577) 
   pls_Cpi = 3.14159265359
   pls_KtoeV = 1.0/11604.52
   

   ! collision diameters for all heavy species
   pls_RSCD(1) = 1.0
   pls_RSCD(2) = 2*1.55
   pls_RSCD(3) = 2*1.52
   pls_RSCD(4) = 1.55
   pls_RSCD(5) = 1.52
   pls_RSCD(6) = 1.52+1.55
   pls_RSCD(7) = 1.0
   pls_RSCD(8) = 2*1.55
   pls_RSCD(9) = 2*1.52
   pls_RSCD(10) = 1.52

   ! molar binary masses
   pls_MHSP(1) = 4.0
   pls_MHSP(2) = 28.0   
   pls_MHSP(3) = 32.0
   pls_MHSP(4) = 14.0
   pls_MHSP(5) = 16.0
   pls_MHSP(6) = 16.0 + 14.0
   pls_MHSP(7) = 4.0
   pls_MHSP(8) = 28.0
   pls_MHSP(9) = 32.0
   pls_MHSP(10)= 16.0

   ! molar mass NaCl + H2O
   pls_MMix = 35.453 + 22.989769 + 18.0

   ! plasma feed rate from jet
   pls_NJET(1) = 0.90*1e26
   pls_NJET(2) = 0.10*0.80*1e26
   pls_NJET(3) = 0.10*0.20*1e26
   pls_NJET(4) = 0.0
   pls_NJET(5) = 0.0
   pls_NJET(6) = 0.0
   pls_NJET(7) = 0.90*1e18
   pls_NJET(8) = 0.10*0.80*1e18
   pls_NJET(9) = 0.10*0.20*1e18
   pls_NJET(10)= 0.01*1e18

   !polynomial coefficients for plasma feed signal
   !done for time and ne without exponential multipliers
   pls_poly_coef(1) = -2.10573490e-12 !deg=21  
   pls_poly_coef(2) =  2.03425802e-10  
   pls_poly_coef(3) = -8.35852479e-09   
   pls_poly_coef(4) =  1.81617147e-07
   pls_poly_coef(5) = -1.92312862e-06  
   pls_poly_coef(6) = -5.59439123e-07   
   pls_poly_coef(7) =  2.56926246e-04  
   pls_poly_coef(8) = -1.97710134e-03
   pls_poly_coef(9) = -2.00464277e-02   
   pls_poly_coef(10)=  4.39665324e-01  
   pls_poly_coef(11)= -1.26430164e+00  
   pls_poly_coef(12)= -4.85148152e+01
   pls_poly_coef(13)=  8.47468391e+02  
   pls_poly_coef(14)= -7.64615720e+03   
   pls_poly_coef(15)=  4.58111902e+04  
   pls_poly_coef(16)= -1.94151558e+05
   pls_poly_coef(17)=  5.89224888e+05  
   pls_poly_coef(18)= -1.25902363e+06   
   pls_poly_coef(19)=  1.80413835e+06  
   pls_poly_coef(20)= -1.55885637e+06
   pls_poly_coef(21)=  6.14457398e+05 !deg=0
   
  if (pls_meshMe .eq. MASTER_PE) then
     write(*,*) 'pls_cfl   =',pls_cfl
     write(*,*) 'pls_sigma =',pls_sigma
     write(*,*) 'pls_dtspec=',pls_dtspec
     write(*,*) 'pls_pct_noise=',pls_pct_noise
  endif

  do lb=1,blockCount


     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)


     plasma_source = nrel
     call Plasma_Feed(plasma_source,rand_noise,                           &
                      solnData(FEED_VAR,:,:,:),solnData(DFUN_VAR,:,:,:),  &
                      blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),         &
                      blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))

      
     solnData(DELE_VAR,:,:,:) = solnData(DELE_VAR,:,:,:) + solnData(FEED_VAR,:,:,:)*pls_dtspec

     do i=0,5
        !feed rate from source
        plasma_source = pls_NJET(i+1)
        call Plasma_Feed(plasma_source,rand_noise,                           &
                         solnData(FEED_VAR,:,:,:),solnData(DFUN_VAR,:,:,:),  &
                         blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),         &
                         blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))

        solnData(DHV0_VAR+i,:,:,:) = solnData(DHV0_VAR+i,:,:,:) + solnData(FEED_VAR,:,:,:)*pls_dtspec
     end do

     do i=0,3
        !feed rate from source
        plasma_source = pls_NJET(i+7)
        call Plasma_Feed(plasma_source,rand_noise,&
                         solnData(FEED_VAR,:,:,:),solnData(DFUN_VAR,:,:,:),&
                         blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                         blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))

        solnData(DHV6_VAR+i,:,:,:) = solnData(DHV6_VAR+i,:,:,:) + solnData(FEED_VAR,:,:,:)*pls_dtspec
     end do
 
     !sum all ions
     solnData(DNIT_VAR,:,:,:) = 0.0
     do i = 0,2
       call Plasma_sumIons(solnData(DHV6_VAR+i,:,:,:),solnData(DNIT_VAR,:,:,:),&
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))
     end do
   
     !find net charge density
     call Plasma_netCharge(solnData(DQNT_VAR,:,:,:),solnData(DNIT_VAR,:,:,:),&
                           solnData(DHV9_VAR,:,:,:),solnData(DELE_VAR,:,:,:),&
                           blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                           blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))
  
     solnData(EPOT_VAR,:,:,:) = 0.0  
     !find induced electric potential
     oldPhi = solnData(EPOT_VAR,:,:,:)
     call Plasma_elPotential(oldPhi,solnData(DFUN_VAR,:,:,:),&
                             solnData(EPOT_VAR,:,:,:),solnData(DQNT_VAR,:,:,:),&
                             del(DIR_X),del(DIR_Y),&
                             blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                             blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS) )

        ! Release block pointers
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  end do

  gcMask = .FALSE.

  gcMask(DELE_VAR)  = .TRUE.
  gcMask(EPOT_VAR)  = .TRUE.

  do i=0,9
        gcMask(DHV0_VAR+i) = .TRUE.
  end do

  call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

  do lb=1,blockCount

     blockID = blockList(lb)

     call Grid_getDeltas(blockID,del)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

     call Plasma_getNorm(solnData(NRMX_VAR,:,:,:), &
                         solnData(NRMY_VAR,:,:,:), &
                         solnData(DFUN_VAR,:,:,:), &
                         blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                         blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                         del(DIR_X),del(DIR_Y))

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

  end do

  gcMask = .FALSE.
  gcMask(NRMX_VAR) = .TRUE.
  gcMask(NRMY_VAR) = .TRUE.

  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

  do lb=1,blockCount

     blockID = blockList(lb)

     call Grid_getDeltas(blockID,del)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

     do j=blkLimits(LOW,JAXIS)-1,blkLimits(HIGH,JAXIS)+1
        do i=blkLimits(LOW,IAXIS)-1,blkLimits(HIGH,IAXIS)+1

           rrc = 0.5 -  solnData(DFUN_VAR,i,j,1)
           rrx = 0.5 - (solnData(DFUN_VAR,i,j,1)+solnData(DFUN_VAR,i-1,j,1))/2.0d0
           rry = 0.5 - (solnData(DFUN_VAR,i,j,1)+solnData(DFUN_VAR,i,j-1,1))/2.0d0

           trc = 0.0
           vrx = 0.0
           vry = 0.0

           do pcount=1,fac
                trc = trc + (pt(pcount)*(rrc**(fac-pcount)));
                vrx = vrx + (pu(pcount)*(rrx**(fac-pcount)));
                vry = vry + (pu(pcount)*(rry**(fac-pcount)));
           end do    

           solnData(TEMP_VAR,i,j,1) = trc*tempDT + Troom

           facexData(VELC_FACE_VAR,i,j,1) = vrx*(solnData(NRMX_VAR,i,j,1)+solnData(NRMX_VAR,i-1,j,1))/2.0d0
           faceyData(VELC_FACE_VAR,i,j,1) = vry*(solnData(NRMY_VAR,i,j,1)+solnData(NRMY_VAR,i,j-1,1))/2.0d0

        end do
     end do

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

  end do

  gcMask = .FALSE.
  gcMask(TEMP_VAR) = .TRUE.
  gcMask(NUNK_VARS+VELC_FACE_VAR) = .TRUE.
  gcMask(NUNK_VARS+1*NFACE_VARS+VELC_FACE_VAR) = .TRUE.

  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
           
  end subroutine Plasma_init
