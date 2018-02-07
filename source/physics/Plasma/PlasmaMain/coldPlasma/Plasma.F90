subroutine Plasma( blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

   use Plasma_data
   use Plasma_interface, only: Plasma_Solve, Plasma_hvDiffCoeff,     &
                               Plasma_elDiffCoeff,Plasma_ColFreq,    &
                               Plasma_sumNeutrals,Plasma_spReactions,&
                               Plasma_spGeneration,Plasma_sumIons,   &
                               Plasma_Feed, Plasma_netCharge,        &
                               Plasma_elPotential

   use Grid_interface, only: Grid_getDeltas, Grid_getBlkIndexLimits,&
                             Grid_getBlkPtr, Grid_releaseBlkPtr,    &
                             Grid_fillGuardCells, Grid_solvePoisson, &
                             GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
                             GRID_PDE_BND_DIRICHLET

   use Driver_data, only: dr_simTime
 
   implicit none
#include "constants.h"
#include "Plasma.h"
#include "Flash.h"   

!#define DEBUG_PLASMA
!#define DEBUG_POISSON

   include "Flash_mpi.h"

   integer, INTENT(INOUT) :: blockCount
   integer, INTENT(IN) :: sweepOrder
   integer, INTENT(INOUT) :: blockList(MAXBLOCKS)
   real,    INTENT(IN) :: timeEndAdv, dt, dtOld

   integer ::  blockID,lb
   real, dimension(MDIM)  :: coord,bsize,del
   real ::  boundBox(2,MDIM)
   integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
   logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
   real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

   
   real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: oldT, oldPhi

   real :: T_res1,T_res,T_resBlock
   real, dimension(10) :: T_resHV, T_resBlockHV
   integer :: ierr, i,j,k
   
   real, dimension(1) :: rand_noise
   real :: plasma_source    !scalar, plasma source rate for species m-3 s-1
   real :: nrel, nrna, nrni !nrh0, nrh1, nrh2, nrh6, nrh7, nrh8, nrh9 
   real :: xcell,ycell,poisfact

   !variables for feed rate
   logical, save :: pls_bullet_flg = .false.
   logical, save :: pls_feed_flg = .false.
   real, save :: pls_bullet_freq = 50e-6
   real, save :: pls_feed_start = 3e-6
   real, save :: pls_feed_end = 12e-6
   real, save :: pls_bullet_timeStamp = 0.0
   real, save :: pls_feed_timeStamp = 0.0
   real, save :: pls_feed_gate = 0.0
   real, save :: pls_freq_gate = 0.0
   real, save :: feed_rate = 0.0
   !end variable for feed rate

   real, parameter :: pi = acos(-1.0)

   integer, dimension(6) :: bc_types
   real, dimension(2,6)  :: bc_values = 0. 

   nrel = 0.99*1e18       ! Electrons
   nrna = 1e26            ! neutral
   nrni = 1e18            ! ions
   plasma_source = 0.0    ! avg val of species in jet 

   T_resBlock      = 0.0
   T_resBlockHV(:) = 0.0

   ! Poisson BCs
   bc_types(1) = GRID_PDE_BND_NEUMANN ! X - Low
   bc_types(2) = GRID_PDE_BND_NEUMANN ! X - High
   bc_types(3) = GRID_PDE_BND_NEUMANN ! Y - Low
   bc_types(4) = GRID_PDE_BND_NEUMANN ! Y - High

   !bc_types(1) = GRID_PDE_BND_DIRICHLET ! X - Low
   !bc_types(2) = GRID_PDE_BND_DIRICHLET ! X - High
   !bc_types(3) = GRID_PDE_BND_DIRICHLET ! Y - Low
   !bc_types(4) = GRID_PDE_BND_DIRICHLET ! Y - High

   do lb = 1,blockCount

     blockID = blockList(lb)

     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)

     call Grid_getBlkCenterCoords(blockId,coord)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

     
     !add all neutral species   
#ifdef DEBUG_PLASMA
     print *,"Going into sumNeurtral"
#endif
     solnData(DNAT_VAR,:,:,:) = 0.0
     do i = 0,5
       call Plasma_sumNeutrals( solnData(DHV0_VAR+i,:,:,:),solnData(DNAT_VAR,:,:,:), &
                                blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                                blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))
     end do 
 
     !add all ions (+ and -)
#ifdef DEBUG_PLASMA
     print *,"Going into sumIons"
#endif
     solnData(DNIT_VAR,:,:,:) = 0.0
     do i = 0,2
       call Plasma_sumIons(solnData(DHV6_VAR+i,:,:,:),solnData(DNIT_VAR,:,:,:),&
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))
     end do

     solnData(DNIT_VAR,:,:,:) = solnData(DNIT_VAR,:,:,:) - solnData(DHV9_VAR,:,:,:)

     !find net charge density in domain
#ifdef DEBUG_PLASMA
     print *,"Going into netCharge"
#endif

     solnData(EPOT_VAR,:,:,:) = 0.0

#ifdef DEBUG_POISSON
     do j=1,blkLimitsGC(HIGH,JAXIS)
        do i=1,blkLimitsGC(HIGH,IAXIS)

           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS)

           ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS)  +  &
                   0.5*del(JAXIS)

           solnData(DQNT_VAR,i,j,1) = -8*pi*pi*cos(2*pi*xcell)*cos(2*pi*ycell)
           !solnData(DQNT_VAR,i,j,1) = -8*pi*pi*sin(2*pi*xcell)*sin(2*pi*ycell)

        end do
     end do
#else
       call Plasma_netCharge(solnData(DQNT_VAR,:,:,:),solnData(DNIT_VAR,:,:,:),&
                             solnData(DHV9_VAR,:,:,:),solnData(DELE_VAR,:,:,:),&
                             blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                             blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))

#endif

     !obtain diffusion coefficient of neutral species
#ifdef DEBUG_PLASMA
     print *,"Going into hvDiffCoeff"
#endif
     do i = 0,5
       call Plasma_hvDiffCoeff(solnData(DFH0_VAR+i,:,:,:), &
                               solnData(PRHV_VAR,:,:,:),solnData(TPHV_VAR,:,:,:), &
                               blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                               blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                               pls_RSCD(i+1), pls_MHSP(i+1)) 
     end do
     
     !obtain collision frequencies
#ifdef DEBUG_PLASMA
     print *,"Going into ColFreq"
#endif
     call Plasma_ColFreq(solnData(FVEI_VAR,:,:,:), solnData(FVEA_VAR,:,:,:), &
                         solnData(DNAT_VAR,:,:,:), solnData(DELE_VAR,:,:,:), &
                         solnData(TPEL_VAR,:,:,:), &
                         blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                         blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS) )


     !obtain diffusion coefficient of electrons (and ions)
#ifdef DEBUG_PLASMA
     print *,"Going into elDiffCoeff"
#endif
     call Plasma_elDiffCoeff(solnData(DFIO_VAR,:,:,:),solnData(DFEL_VAR,:,:,:),&
                             solnData(TPEL_VAR,:,:,:),solnData(TPHV_VAR,:,:,:),&
                             solnData(FVEA_VAR,:,:,:),solnData(FVEI_VAR,:,:,:),&
                             blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                             blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS) )
      
     !obtain chemical reaction rates for all species
#ifdef DEBUG_PLASMA
     print *,"Going into spReactions"
#endif
     call Plasma_spReactions(solnData(RSP0_VAR,:,:,:), solnData(RSP1_VAR,:,:,:),&
                             solnData(RSP2_VAR,:,:,:), solnData(RSP3_VAR,:,:,:),&
                             solnData(RSP4_VAR,:,:,:), solnData(RSP5_VAR,:,:,:),&
                             solnData(RSP6_VAR,:,:,:), solnData(RSP7_VAR,:,:,:),&
                             solnData(RSP8_VAR,:,:,:), solnData(RSP9_VAR,:,:,:),&
                             solnData(RSP10_VAR,:,:,:), solnData(RSP11_VAR,:,:,:),&
                             solnData(RSP12_VAR,:,:,:), solnData(RSP13_VAR,:,:,:),&
                             solnData(RSP14_VAR,:,:,:), solnData(RSP15_VAR,:,:,:),&
                             solnData(RSP16_VAR,:,:,:),&
                             solnData(TPHV_VAR,:,:,:),  solnData(TPEL_VAR,:,:,:),&
                             blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                             blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS) )

     !obtain species generation rates using chemical reations
#ifdef DEBUG_PLASMA
     print *,"Going into spGeneration"
#endif
     call Plasma_spGeneration(solnData(DHV0_VAR,:,:,:), solnData(DHV1_VAR,:,:,:),& 
                              solnData(DHV2_VAR,:,:,:), solnData(DHV3_VAR,:,:,:),&
                              solnData(DHV4_VAR,:,:,:), solnData(DHV5_VAR,:,:,:),&
                              solnData(DHV6_VAR,:,:,:), solnData(DHV7_VAR,:,:,:),&
                              solnData(DHV8_VAR,:,:,:), solnData(DHV9_VAR,:,:,:),&
                              solnData(DELE_VAR,:,:,:),                          &
                              solnData(RSP0_VAR,:,:,:), solnData(RSP1_VAR,:,:,:),&
                              solnData(RSP2_VAR,:,:,:), solnData(RSP3_VAR,:,:,:),&
                              solnData(RSP4_VAR,:,:,:), solnData(RSP5_VAR,:,:,:),&
                              solnData(RSP6_VAR,:,:,:), solnData(RSP7_VAR,:,:,:),&
                              solnData(RSP8_VAR,:,:,:), solnData(RSP9_VAR,:,:,:),&
                              solnData(RSP10_VAR,:,:,:),solnData(RSP11_VAR,:,:,:),&
                              solnData(RSP12_VAR,:,:,:),solnData(RSP13_VAR,:,:,:),&
                              solnData(RSP14_VAR,:,:,:), solnData(RSP15_VAR,:,:,:),&
                              solnData(RSP16_VAR,:,:,:),&
                              solnData(GNH0_VAR,:,:,:), solnData(GNH1_VAR,:,:,:),& 
                              solnData(GNH2_VAR,:,:,:), solnData(GNH3_VAR,:,:,:),&
                              solnData(GNH4_VAR,:,:,:), solnData(GNH5_VAR,:,:,:),&
                              solnData(GNH6_VAR,:,:,:), solnData(GNH7_VAR,:,:,:),&
                              solnData(GNH8_VAR,:,:,:), solnData(GNH9_VAR,:,:,:),&
                              solnData(GNE_VAR,:,:,:),  solnData(GNEBZ_VAR,:,:,:),&
                              solnData(GNERT_VAR,:,:,:),                          &
                              blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                              blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS) )     
 
     plasma_source = 0.0
     !obtain number density of electrons
#ifdef DEBUG_PLASMA
     print *,"Going into feed rate of charged particles"
#endif

     !*****************************************************************************!     
     !*****************************************************************************!
     !plasma jet feed rate for electrons and ions

     !Check if bullet frequency criteria met
     pls_freq_gate = dr_simTime - pls_bullet_timeStamp 
    
     if( (pls_bullet_flg .eqv. .false.) .and. (pls_freq_gate.ge.pls_bullet_freq) ) then
        pls_feed_timeStamp = dr_simTime
        pls_bullet_flg = .true.
     else
        pls_feed_timeStamp = pls_feed_timeStamp
        pls_bullet_flg = .false.
     end if

     !Check if feed rate criteria met
     !pls_feed_gate = dr_simTime - pls_feed_timeStamp

     if(pls_bullet_flg .eqv. .true.) then
 
        pls_feed_gate = dr_simTime - pls_feed_timeStamp

        if( (pls_feed_flg .eqv. .false.) .and. (pls_feed_gate.ge.pls_feed_start) ) then
           pls_feed_flg = .true.
        else
           pls_feed_flg = .false.
        end if

        !Check if feed rate can be computed
        if (pls_feed_flg .eqv. .true.) then

           if((pls_feed_gate.ge.pls_feed_start).and.&
              (pls_feed_gate.le.pls_feed_end)) then

              plasma_source = 0.0
              ! nom nom      
              do j=1,21
                 plasma_source = plasma_source + & 
                                 1e18*pls_poly_coef(j)*(1e6*pls_feed_gate**(21-j))
              end do

           else
              plasma_source = 0.0
         
           end if    
        
        end if

        !Initialize values for next bullet
        if (pls_feed_gate.gt.pls_feed_end) then
           plasma_source = 0.0
           pls_bullet_flg = .false.
           pls_feed_flg = .false.
           pls_bullet_timeStamp = dr_simTime
        end if

        !end if    

     end if    

     !*****************************************************************************!
     !*****************************************************************************! 

     call Plasma_Feed(plasma_source,rand_noise,                           &
                      solnData(FEED_VAR,:,:,:),solnData(DFUN_VAR,:,:,:),  &
                      blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),         &
                      blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))

     !number density of electrons
#ifdef DEBUG_PLASMA
     print *,"Going into Plasma solve electrons"
#endif

     !reference density value
     oldT = solnData(DELE_VAR,:,:,:)
     !solve for new density
     call Plasma_Solve(solnData(DELE_VAR,:,:,:), solnData(GNE_VAR,:,:,:), & 
                       oldT, solnData(FEED_VAR,:,:,:),                    &
                       solnData(DFUN_VAR,:,:,:), solnData(DFEL_VAR,:,:,:),&
                       dt,del(DIR_X),del(DIR_Y),                          &
                       blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),        &
                       blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),T_res1)
     T_resBlock = T_resBlock + T_res1
      
     !number density of ions
#ifdef DEBUG_PLASMA
     print *,"Going into Plasma solve ions"
#endif
     do i=0,3
        !reference density values
        oldT = solnData(DHV6_VAR+i,:,:,:)
        !solve for new density
        call Plasma_Solve(solnData(DHV6_VAR+i,:,:,:),solnData(GNH6_VAR+i,:,:,:),& 
                          oldT, solnData(FEED_VAR,:,:,:),                       &
                          solnData(DFUN_VAR,:,:,:),  solnData(DFEL_VAR,:,:,:),  &
                          dt,del(DIR_X),del(DIR_Y),                             &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),           &
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),T_res1)
        T_resBlockHV(i+7) = T_resBlockHV(i+7) + T_res1
     end do

     !obtain number density of heavy species
#ifdef DEBUG_PLASMA
     print *,"Going into Plasma solve heavy"
#endif
     do i=0,5
        !feed rate from source, added at every time step
        plasma_source = pls_NJET(i+1)
        call Plasma_Feed(plasma_source,rand_noise,                           &
                         solnData(FEED_VAR,:,:,:),solnData(DFUN_VAR,:,:,:),  &
                         blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),         &
                         blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))
        !reference density values
        oldT = solnData(DHV0_VAR+i,:,:,:)
        !sove for new density
        call Plasma_Solve(solnData(DHV0_VAR+i,:,:,:), solnData(GNH0_VAR+i,:,:,:),&
                          oldT, solnData(FEED_VAR,:,:,:),                        &
                          solnData(DFUN_VAR,:,:,:), solnData(DFH0_VAR+i,:,:,:),&
                          dt,del(DIR_X),del(DIR_Y),                              &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),            &
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),T_res1)
        T_resBlockHV(i+1) = T_resBlockHV(i+1) + T_res1
     end do
     
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)


   end do

   if(blockCount .gt. 0) T_resBlock   = T_resBlock/blockCount
   if(blockCount .gt. 0) T_resBlockHV = T_resBlockHV/blockCount

   ! Collect residuals from other processes
   call MPI_Allreduce(T_resBlock, T_res, 1, FLASH_REAL,&
                      MPI_SUM, MPI_COMM_WORLD, ierr)
   
   call MPI_Allreduce(T_resBlockHV, T_resHV, 10, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   T_res   = T_res/pls_meshNumProcs
   T_resHV = T_resHV/pls_meshNumProcs

   if(pls_meshMe .eq. MASTER_PE) then
      
        print *,"Ne residual","",":",T_res

        do i=1,10
          print *,"Nh residual",i,":",T_resHV(i)
        end do

   endif

  gcMask = .FALSE.

  gcMask(DELE_VAR)  = .TRUE.
  gcMask(DFEL_VAR)  = .TRUE.
  gcMask(FVEA_VAR)  = .TRUE.
  gcMask(FVEI_VAR)  = .TRUE.
  gcMask(DNAT_VAR)  = .TRUE.
  gcMask(DNIT_VAR)  = .TRUE. 
  gcMask(TPHV_VAR)  = .TRUE.
  gcMask(TPEL_VAR)  = .TRUE.
  gcMask(GNE_VAR)   = .TRUE.
  gcMask(GNEBZ_VAR) = .TRUE.
  gcMask(GNERT_VAR) = .TRUE.

  gcMask(DQNT_VAR) = .TRUE.

  do i=0,9
        gcMask(DHV0_VAR+i) = .TRUE.
        gcMask(DFH0_VAR+i) = .TRUE. 
        gcMask(GNH0_VAR+i) = .TRUE. 
  end do

  do i=0,16
        gcMask(RSP0_VAR+i) = .TRUE.
  end do

  call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)


  poisfact = 1.0
  call Grid_solvePoisson (EPOT_VAR, DQNT_VAR, bc_types, bc_values, poisfact)

  gcMask = .FALSE.
  gcMask(EPOT_VAR) = .TRUE.

  call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
  

end subroutine Plasma
