!!****if* source/Simulation/SimulationMain/Nuc2Grid/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  call Driver_evolveFlash()
!!
!! DESCRIPTION
!!
!!  Pseudo evolution for Phoenix input data converter.
!!
!!
!!***

!#define DEBUG_DRIVER
#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

!!REORDER(4):solnData

subroutine Driver_evolveFlash()

  use Driver_data,         ONLY : dr_globalMe, dr_nbegin,       &
                                  dr_nend, dr_dt, dr_wallClockTimeLimit, &
                                  dr_tmax, dr_simTime, dr_redshift,      &
                                  dr_nstep, dr_dtOld, dr_dtNew,          &
                                  dr_simGeneration,                      &
                                  dr_restart, dr_elapsedWCTime,          &
                                  dr_redshiftOld, dr_useRedshift,        &
                                  dr_redshiftfinal
  use Simulation_data, ONLY : sim_doWAvg,sim_doGP,sim_doWeight, &
       sim_doConvolve,sim_doInterpExtrap, &
       sim_doLowerBounds,sim_doEos,&
       sim_smlrho,sim_smallE,sim_smallT, sim_abundanceFixupMaxDens, &
       sim_doFixupAbundances, sim_unkCellWeight, sim_useTrajValues
  use sim_outputGridData,  ONLY:  ndimOg, geometryOg
  use Grid_data,           ONLY:  gr_vartypes
  use Driver_interface,    ONLY : Driver_sourceTerms, Driver_computeDt, &
                                  Driver_getElapsedWCTime, Driver_abortFlash
  use Logfile_interface,   ONLY : Logfile_stamp, Logfile_close
  use Timers_interface,    ONLY : Timers_start, Timers_stop, &
                                  Timers_getSummary
  use Particles_interface, ONLY : Particles_advance
  use Grid_interface,      ONLY : Grid_getLocalNumBlks, &
                                  Grid_getListOfBlocks, &
                                  Grid_updateRefinement,&
                                  Grid_getBlkIndexLimits,&
                                  Grid_getBlkPtr,&
                                  Grid_releaseBlkPtr,&
                                  Grid_fillGuardCells, Grid_getCellCoords,&
                                  Grid_dump
  use Eos_interface,       ONLY : Eos_wrapped
  use Hydro_interface,     ONLY : Hydro
  use Gravity_interface,   ONLY : Gravity_potentialListOfBlocks
  use IO_interface,        ONLY : IO_output,IO_outputFinal
  use Cosmology_interface, ONLY:  Cosmology_redshiftHydro, &
                                  Cosmology_solveFriedmannEqn, &
                                  Cosmology_getRedshift


  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer   :: localNumBlocks

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)
  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)
  integer       :: bcTypes(6)
  real          :: bcValues(2,6) = 0.
  integer       :: ivar, iKern, lb, i,j,k

  integer :: ii  ! test, by Min

  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(4,2) :: strBuff
  character(len=15) :: numToStr

  logical :: gridChanged
  logical :: endRun
!  logical :: skip

  real    :: dvol
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  integer :: point(MDIM), error
  integer,parameter :: nGlobalSum = 1
  real :: gsum(nGlobalSum) !Global summed quantities
  real :: lsum(nGlobalSum) !Global summed quantities
  real :: renormFact(nGlobalSum)
  real :: sumZiYi, weightFromNumP
!  integer :: species2Z(SPECIES_BEGIN:SPECIES_END)
  real,allocatable,dimension(:) :: xInterp, yInterp
  integer,allocatable,dimension(:) :: indInterp


!  data species2Z /13,18, 5,4, 6,20,17,27,24,29, 9,26, 31,32, 1,2, 19, 3, 12,25, 7,11,10, 0, 28, 8, 15, 16,21,14, 22, 23, 30/
  endRun = .false.

  call Logfile_stamp('Entering evolution routine' , '[Driver_evolveFlash]')

  iKern = GAUS_VAR

  call Grid_getLocalNumBlks(localNumBlocks)
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  ! Sum quantities over all locally held leaf-node blocks.
  gsum  = 0.
  lsum = 0.
  do lb = 1, blockCount
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)

     ! Sum contributions from the indicated blkLimits of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k

              !! Get the cell volume for a single cell
!!$              call Grid_getSingleCellVol(blockList(lb), EXTERIOR, point, dvol)
              dvol = 1.0

              ! smeared particle count
              lsum(1) = lsum(1) + solnData(iKern,i,j,k)*dvol 
           enddo
        enddo
     enddo
     call Grid_releaseBlkPtr(blockList(lb), solnData)

  enddo



  ! Now every PE sums the local contributions from all of
  ! the processors

  call MPI_ALLREDUCE (lsum, gsum, nGlobalSum, FLASH_REAL, MPI_SUM, & 
       &                MPI_COMM_WORLD, error)

  if (gsum(1) .NE. 0) then
     renormFact(1) = 1.0 / gsum(1)
     do lb = 1, blockCount
        !get the index limits of the block
        call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

        ! get a pointer to the current block of data
        call Grid_getBlkPtr(blockList(lb), solnData)

        ! Sum contributions from the indicated blkLimits of cells.
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 solnData(iKern,i,j,k) = solnData(iKern,i,j,k) * renormFact(1)
              enddo
           enddo
        enddo
        call Grid_releaseBlkPtr(blockList(lb), solnData)
     end do
  end if
  

  call Logfile_stamp('Entering evolution loop' , '[Driver_evolveFlash]')
  call Timers_start("evolution")

  dr_simGeneration = 0

  do dr_nstep = dr_nBegin, dr_nend
     print*, 'do dr_nstep = dr_nBegin, dr_nend ;  ', dr_nstep, dr_nBegin, dr_nend
     
     !!Step forward in time. See bottom of loop for time step calculation.
     call Grid_getLocalNumBlks(localNumBlocks)
     call Grid_getListOfBlocks(LEAF,blockList,blockCount)

     call doLog


     dr_simGeneration = dr_simGeneration + 1 

     !--------------------------------------------------------------------
     !- Start Physics Sequence
     !----

     select case (dr_nstep - dr_nBegin + 1)
     case(1)
!!$        print*,'Converting TEMP_VAR from [10^9 K] to [K]'
!!$        do lb = 1, blockCount
!!$           call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)
!!$           call Grid_getBlkPtr(blockList(lb), solnData)
!!$           do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
!!$              do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
!!$                 do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
!!$                    solnData(TEMP_VAR,i,j,k) = solnData(TEMP_VAR,i,j,k) * 1.0e9
!!$                 end do
!!$              end do
!!$           end do
!!$           call Grid_releaseBlkPtr(blockList(lb), solnData)
!!$        end do

!!$        do lb = 1, blockCount
!!$           call Grid_dump((/DENS_VAR/),1,blockList(lb),gcell=.FALSE.)
!!$        end do

        if (sim_doConvolve) then
           print*,'Convolving UNK var #',NUMP_VAR,' (NUMP_VAR)'
           call Grid_convolve(NUMP_VAR,NUMP_VAR,iKern,bcTypes,bcValues,1.0)
           do lb = 1, blockCount
              !get the index limits of the block
              call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

              ! get a pointer to the current block of data
              call Grid_getBlkPtr(blockList(lb), solnData)

              ! Sum contributions from the indicated blkLimits of cells.
              do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                 do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                    do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                       solnData(NUMP_VAR,i,j,k) = max(solnData(NUMP_VAR,i,j,k),1.0e-6)
                    enddo
                 enddo
              enddo
              call Grid_releaseBlkPtr(blockList(lb), solnData)

           enddo
        else
           print*,'NOT convolving UNK var #',NUMP_VAR,' (NUMP_VAR)'
        end if

        do lb = 1, blockCount
           !get the index limits of the block
           call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

           ! get a pointer to the current block of data
           call Grid_getBlkPtr(blockList(lb), solnData)

           ! Sum contributions from the indicated blkLimits of cells.
           do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
              do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                 do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                    solnData(NUP1_VAR,i,j,k) =  solnData(NUMP_VAR,i,j,k)
                 enddo
              enddo
           enddo
           call Grid_releaseBlkPtr(blockList(lb), solnData)

        enddo

     case(2)
        if (sim_doConvolve) then
           do ivar = 1,NUNK_VARS
              if (ivar  .NE. iKern .AND. &
                   ivar .NE. PDEN_VAR .AND. &
                   ivar .NE. NUP0_VAR .AND. &
                   ivar .NE. NUP1_VAR .AND. &
                   ivar .NE. NUMP_VAR .AND. &
                   ivar .NE. DENS_VAR) then
                 print*,'Convolving UNK var #',ivar
#define SMEARING_CONSERVATIVE .true.
                 if (SMEARING_CONSERVATIVE .AND. (gr_vartypes(ivar).eq.VARTYPE_PER_MASS)) then
                    print*,'Rescaling UNK var #',ivar,' with DENS_VAR'
                    do lb = 1, blockCount
                       !get the index limits of the block
                       call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

                       ! get a pointer to the current block of data
                       call Grid_getBlkPtr(blockList(lb), solnData)

                       ! Sum contributions from the indicated blkLimits of cells.
                       do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                          do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                             do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                                if (sim_useTrajValues(DENS_VAR)) then
                                   weightFromNumP = solnData(NUP0_VAR,i,j,k) + sim_unkCellWeight(DENS_VAR)
                                else
                                   weightFromNumP = 1.0
                                end if
                                if (weightFromNumP.NE.0.0) then
                                   solnData(ivar,i,j,k) = solnData(ivar,i,j,k) * &
                                        solnData(DENS_VAR,i,j,k) / weightFromNumP
                                else
                                   solnData(ivar,i,j,k) = 0.0
                                end if
                             enddo
                          enddo
                       enddo
                       call Grid_releaseBlkPtr(blockList(lb), solnData)
                    enddo
                 end if

                 call Grid_convolve(ivar,ivar,iKern,bcTypes,bcValues,1.0)
              end if
           end do
           ivar = DENS_VAR
           print*,'Convolving UNK var #',ivar,' (DENS_VAR)'
           call Grid_convolve(ivar,ivar,iKern,bcTypes,bcValues,1.0)
           if(sim_useTrajValues(ivar)) then
              print*,'Rescaling UNK var #',ivar,' (DENS_VAR) with NUMP_VAR'
              do lb = 1, blockCount
                 !get the index limits of the block
                 call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

                 ! get a pointer to the current block of data
                 call Grid_getBlkPtr(blockList(lb), solnData)

                 ! Sum contributions from the indicated blkLimits of cells.
                 do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                    do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                       do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                          weightFromNumP = solnData(NUMP_VAR,i,j,k) + sim_unkCellWeight(ivar)
                          solnData(ivar,i,j,k) = solnData(ivar,i,j,k) / weightFromNumP
                       enddo
                    enddo
                 enddo
                 call Grid_releaseBlkPtr(blockList(lb), solnData)
              enddo
           else
              print*,'NOT Rescaling UNK var #',ivar,' (DENS_VAR) with NUMP_VAR, since not using trajectory values'
           end if
           do ivar = 1,NUNK_VARS
              if (ivar  .NE. iKern .AND. &
                   ivar .NE. PDEN_VAR .AND. &
                   ivar .NE. NUP0_VAR .AND. &
                   ivar .NE. NUP1_VAR .AND. &
                   ivar .NE. NUMP_VAR .AND. &
                   ivar .NE. DENS_VAR) then
                 if(sim_useTrajValues(ivar) .OR. (SMEARING_CONSERVATIVE .AND. (gr_vartypes(ivar).eq.VARTYPE_PER_MASS))) then
                    if (.NOT. sim_useTrajValues(ivar)) then
                       print*,'Rescaling UNK var #',ivar,' with DENS_VAR'
                    else if (SMEARING_CONSERVATIVE .AND. (gr_vartypes(ivar).eq.VARTYPE_PER_MASS)) then
                       print*,'Rescaling UNK var #',ivar,' with NUMP_VAR and DENS_VAR'
                    else
                       print*,'Rescaling UNK var #',ivar,' with NUMP_VAR'
                    end if
                    do lb = 1, blockCount
                       !get the index limits of the block
                       call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

                       ! get a pointer to the current block of data
                       call Grid_getBlkPtr(blockList(lb), solnData)

                       ! Sum contributions from the indicated blkLimits of cells.
                       do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                          do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                             do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                                if (SMEARING_CONSERVATIVE .AND. (gr_vartypes(ivar).eq.VARTYPE_PER_MASS)) then
                                   solnData(ivar,i,j,k) = solnData(ivar,i,j,k) / solnData(DENS_VAR,i,j,k)
                                end if
                                if (sim_useTrajValues(ivar)) then
                                   weightFromNumP = solnData(NUMP_VAR,i,j,k) + sim_unkCellWeight(ivar)
                                   solnData(ivar,i,j,k) = solnData(ivar,i,j,k) / weightFromNumP
                                end if
                             enddo
                          enddo
                       enddo
                       call Grid_releaseBlkPtr(blockList(lb), solnData)

                    enddo
                 end if
              end if
           end do


        else
           print*,'NOT convolving other UNK vars #'
           do ivar = 1,NUNK_VARS
              if (ivar  .NE. iKern .AND. &
                   ivar .NE. PDEN_VAR .AND. &
                   ivar .NE. NUP0_VAR .AND. &
                   ivar .NE. NUP1_VAR .AND. &
                   ivar .NE. NUMP_VAR) then
                 if (sim_useTrajValues(ivar)) then
                    print*,'Rescaling UNK var #',ivar,' with NUMP_VAR'
                    do lb = 1, blockCount
                       !get the index limits of the block
                       call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

                       ! get a pointer to the current block of data
                       call Grid_getBlkPtr(blockList(lb), solnData)

                       ! Sum contributions from the indicated blkLimits of cells.
                       do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                          do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                             do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                                weightFromNumP = solnData(NUMP_VAR,i,j,k) + sim_unkCellWeight(ivar)
                                if (weightFromNumP .NE. 0.0) then
                                   solnData(ivar,i,j,k) = solnData(ivar,i,j,k) / weightFromNumP
#ifdef DEBUG_DRIVER
                                else
                                   print*,'Not scaling',ivar,': NUMP is 0.0 at i/j=',i,j
#endif
                                end if
                             enddo
                          enddo
                       enddo
                       call Grid_releaseBlkPtr(blockList(lb), solnData)

                    enddo
                 end if
              end if
           end do
        end if
     case(3)
        if (sim_doWeight) then
          ! Fill in empty cells by weighting nearby cells   
          print*,'calling Grid_fillGuardCells...'
          call Grid_fillGuardCells(CENTER,ALLDIR,selectBlockType=LEAF)
          print*,'Applying weighting for any leaf node UNK cells that got no particles - 1D,2D,3D'

          do lb = 1, blockCount
             call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)
             call Grid_getBlkPtr(blockList(lb), solnData)

             ! do weighting 
             call doWeight(blockList(lb))

             call Grid_releaseBlkPtr(blockList(lb), solnData)
          end do ! lb
          print*,'finished doWeight' 
        else
          print*,'NOT applying weighting for UNK cells that got no particles'
        end if


        if (sim_doInterpExtrap) then
        ! Fill in missing cells by interpolation, or extrapolation where necessary.
           print*,'calling Grid_fillGuardCells...'
           ! currently only interpolates/extrapolates in X direction...
           call Grid_fillGuardCells(CENTER,IAXIS,selectBlockType=LEAF)
           print*,'Applying interpolation / extrapolation for any leaf node UNK cells that got no particles - X direction'
           do lb = 1, blockCount
              call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)
              allocate(xInterp(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))
              allocate(yInterp(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))
              allocate(indInterp(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))
              call Grid_getBlkPtr(blockList(lb), solnData)
              do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                 do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                    ! currently only interpolates/extrapolates in X direction...
                    call doInterpExtrap(blockList(lb),IAXIS)
                 end do
              end do
              call Grid_releaseBlkPtr(blockList(lb), solnData)
              deallocate(xInterp)
              deallocate(yInterp)
              deallocate(indInterp)
           end do
           if (NDIM > 1) then
              print*,'calling Grid_fillGuardCells...'
              call Grid_fillGuardCells(CENTER,JAXIS,selectBlockType=LEAF)
              print*,'Applying interpolation / extrapolation for any leaf node UNK cells that got no particles - now Y direction'
              do lb = 1, blockCount
                 call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)
                 allocate(xInterp(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)))
                 allocate(yInterp(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)))
                 allocate(indInterp(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)))
                 call Grid_getBlkPtr(blockList(lb), solnData)
                 do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                    do j = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                       ! Now only interpolates/extrapolates in Y direction...
                       call doInterpExtrap(blockList(lb),JAXIS)
                    end do
                 end do
                 call Grid_releaseBlkPtr(blockList(lb), solnData)
                 deallocate(xInterp)
                 deallocate(yInterp)
                 deallocate(indInterp)
              end do
           end if
           print*,'finished doInterExter'
        else
           print*,'NOT applying interpolation / extrapolation for UNK cells that got no particles'
        end if

     case(4)
        if (sim_doLowerBounds) then
        ! Apply lower bounds on values
           print*,'Applying lower limits to dens,temp,eint'
           do lb = 1, blockCount
              call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)
              call Grid_getBlkPtr(blockList(lb), solnData)
              do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                 do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                    do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                       solnData(TEMP_VAR,i,j,k) = max(solnData(TEMP_VAR,i,j,k),sim_smallT)
                       solnData(EINT_VAR,i,j,k) = max(solnData(EINT_VAR,i,j,k),sim_smallE)
#ifdef ENER_VAR
                       solnData(ENER_VAR,i,j,k) = max(solnData(ENER_VAR,i,j,k),sim_smallE)
#endif
                       solnData(DENS_VAR,i,j,k) = max(solnData(DENS_VAR,i,j,k),sim_smlrho)
                    end do
                 end do
              end do
              call Grid_releaseBlkPtr(blockList(lb), solnData)
           end do
           print*, 'finished applying lower limits'
        else
           print*,'NOT applying lower limits to dens,temp,eint'
        end if


        if (.TRUE.) then
        ! Computing Sum X_i
           print*,'Computing Sum X_i'
           do lb = 1, blockCount
              call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)
              call Grid_getBlkPtr(blockList(lb), solnData)
              do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                 do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                    do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                       sumZiYi = 0.0
                       solnData(SUMX_MSCALAR,i,j,k) = sum(solnData(SPECIES_BEGIN:SPECIES_END,i,j,k))

                       ! check zeros in SolnData -by Min
                       !if (solnData(SUMX_MSCALAR,i,j,k)==0) then
                         !print*, 'SolnData(SUMX_MSCALAR)=0',solnData(SUMX_MSCALAR,i,j,k)
                         !print*, 'soln sum', sum(solnData(SPECIES_BEGIN:SPECIES_END,i,j,k))
                         !do ii =  SPECIES_BEGIN, SPECIES_END
                           !print*, 'soln elem', ii, solnData(ii,i,j,k)
                         !end do
                       !end if
                       ! end of check

                    end do
                 end do
              end do
              call Grid_releaseBlkPtr(blockList(lb), solnData)
           end do
        else
           print*,'NOT computing Sum X_i'
        end if

        if (sim_doFixupAbundances) then
        ! Do fixum where density is below a threshold value
           print*,'Fixing up abundances by renorming sum to 1.'
           do lb = 1, blockCount
              call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)
              call Grid_getBlkPtr(blockList(lb), solnData)
              do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                 do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                    do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                       if (solnData(DENS_VAR,i,j,k) < sim_abundanceFixupMaxDens) then

                          ! check denominator = 0 or not -by Min 
                          if (solnData(SUMX_MSCALAR,i,j,k) <= 1.0e-8) then
                          ! print*, 'solnData(SUMX_MSCALAR,i,j,k) <= 1.0e-8)', solnData(SUMX_MSCALAR,i,j,k)
                           solnData(SUMX_MSCALAR,i,j,k) = 1.0e-8
                          end if
                          ! end of check

                          solnData(SPECIES_BEGIN:SPECIES_END,i,j,k) = &
                               solnData(SPECIES_BEGIN:SPECIES_END,i,j,k) / solnData(SUMX_MSCALAR,i,j,k)
                          solnData(NI56_MSCALAR,i,j,k) = solnData(NI56_MSCALAR,i,j,k) / solnData(SUMX_MSCALAR,i,j,k)
                       end if
                    end do
                 end do
              end do
              call Grid_releaseBlkPtr(blockList(lb), solnData)
           end do
        else
           print*,'NOT fixing up abundances'
        end if
        

   !  case(5)
   !     if (sim_doEos) then
   !        ! Compute pressure
   !        ! print*,'Applying Eos_wrapped, mode',MODE_DENS_TEMP
   !        ! do lb = 1, blockCount
   !        !   call Eos_wrapped(MODE_DENS_TEMP,blkLimits, blockList(lb))
   !        !end do
   !     else
   !        print*,'NOT calling Eos_wrapped'
   !     end if

   !  case(6)
   !     if (ndimOg == 2 .AND. geometryOg == CYLINDRICAL) then
           !print*,'Filling VELZ_VAR from VELY_VAR',MODE_DENS_TEMP
           !do lb = 1, blockCount
           !   call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)
           !   call Grid_getBlkPtr(blockList(lb), solnData)
           !   do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           !      do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           !         do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
           !            solnData(VELZ_VAR,i,j,k) = solnData(VELY_VAR,i,j,k)
           !            solnData(VELY_VAR,i,j,k) = 0.0
           !         end do
           !      end do
           !   end do
           !   call Grid_releaseBlkPtr(blockList(lb), solnData)
           !end do
    !    end if
    !    ! Map to output grid
        call sim_setupOutputGrid(dr_globalMe)
        call sim_mapUnkVarsToOutputGrid
        call sim_shareOutputGrid(dr_globalMe)
    ! case(7)
     case(5)
        if (dr_globalMe == MASTER_PE) then
           ! write the output grid
           call sim_writeOutputGrid
        end if
     end select


     !output a plotfile before the grid changes
     write(*,'(A25,2i5)') 'start to output plt file', dr_nstep, dr_nstep-dr_nBegin+1

     call Timers_start("IO_output")
     call IO_output(dr_simTime, &
          dr_dt, dr_nstep+1, dr_nbegin, endRun, PLOTFILE_AND_PARTICLEFILE)
     call Timers_stop("IO_output")

     write(*,'(A25,2i5)') 'finished outputing plt file', dr_nstep, dr_nstep-dr_nBegin+1

!!$     call Timers_start("Grid_updateRefinement")
!!$     call Grid_updateRefinement(dr_nstep, dr_simTime, gridChanged)
!!$     call Timers_stop("Grid_updateRefinement")


     write(*,'(A25,2i5)') 'start to output chk file', dr_nstep, dr_nstep-dr_nBegin+1

     call Timers_start("IO_output")
     call IO_output(dr_simTime,dr_dt,dr_nstep+1,dr_nbegin,endRun,&
          CHECKPOINT_FILE_ONLY)
     call Timers_stop("IO_output")

     write(*,'(A25,2i5)') 'finished outputing chk file', dr_nstep, dr_nstep-dr_nBegin+1


  enddo
  !The value of dr_nstep after the loop is (dr_nend + 1) if the loop iterated for
  !the maximum number of times.  However, we need to retain the value that
  !dr_nstep had during the last loop iteration, otherwise the number for nstep
  !that will be stored in a final checkpoint file will be wrong.
  dr_nstep = min(dr_nstep,dr_nend)


  call Timers_stop("evolution")
  call Logfile_stamp('Exiting evolution loop' , '[Driver_evolveFlash]')
  if(.NOT.endRun) call IO_outputFinal()

  call Timers_getSummary(max(0,dr_nstep-dr_nbegin+1))


  call Logfile_stamp("FLASH run complete.", "LOGFILE_END")
  call Logfile_close()

  return

contains
  subroutine doLog
    if (dr_globalMe == MASTER_PE) then

       write (numToStr(1:), '(I10)') dr_nstep
       write (strBuff(1,1), "(A)") "n"
       write (strBuff(1,2), "(A)") trim(adjustl(numToStr))

       write (numToStr(1:), "(1PE12.6)") dr_simTime
       write (strBuff(2,1), "(A)") "t"
       write (strBuff(2,2), "(A)") trim(adjustl(numToStr))

       if (.not. dr_useRedshift) then

          write (numToStr(1:), "(1PE12.6)") dr_dt
          write (strBuff(3,1), "(A)") "dt"
          write (strBuff(3,2), "(A)") trim(adjustl(NumToStr))

          call Logfile_stamp(strBuff(1:3,:), 3, 2, "step")

       else

          write (numToStr(1:), "(F8.3)") dr_redshift
          write (strBuff(3,1), "(A)") "z"
          write (strBuff(3,2), "(A)") trim(adjustl(NumToStr))

          write (numToStr(1:), "(1PE12.6)") dr_dt
          write (strBuff(4,1), "(A)") "dt"
          write (strBuff(4,2), "(A)") trim(adjustl(NumToStr))

          call Logfile_stamp(strBuff, 4, 2, "step")

       endif

    end if
  end subroutine doLog

  subroutine doInterpExtrap(blockId,axis)
    integer,intent(in) :: blockId
    integer,intent(in) :: axis
    integer :: i, iInterp, nInterp, ih, holeCount, iLoc
    integer,allocatable :: holeInd(:)
    real,allocatable :: xHole(:)
    integer :: sizeX
    real :: newVal,oldVal, dummydy

    iInterp = 0
    holeCount = 0
    if (axis==IAXIS) then
       ! currently first interpolates/extrapolates in X direction...
       do i = blkLimitsGC(LOW,IAXIS)+1, blkLimitsGC(HIGH,IAXIS)-1
          if (solnData(NUMP_VAR,i,j,k) == 0.0) then
             holeCount = holeCount + 1
          end if
       end do
       if (holeCount==0) return
       allocate(holeInd(holeCount))
       allocate(xHole(holeCount))
       holeCount = 0
       do i = blkLimitsGC(LOW,IAXIS)+1, blkLimitsGC(HIGH,IAXIS)-1
          if (solnData(NUMP_VAR,i,j,k) == 0.0) then
             holeCount = holeCount + 1
             holeInd(holeCount) = i
          end if
       end do
       indInterp(:) = 0
       do ih = 1,holeCount
          indInterp(holeInd(ih)) = ih
       end do

       sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
       call Grid_getCellCoords(IAXIS, blockId, CENTER, .TRUE., xInterp, sizeX)

       iInterp = 0; ih = 0
       do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
          if (indInterp(i) == 0) then
             iInterp = iInterp + 1
             if (iInterp< i) xInterp(iInterp) = xInterp(i)
          else  !  (indInterp(i) > 0)
             ih = ih + 1
             xHole(ih) = xInterp(i)
          end if
       end do
       nInterp = iInterp
       if (nInterp + holeCount .NE. sizeX) then
          print*,"doInterpExtrap: block index confusion; nInterp,holeCount,sizeX=",nInterp,holeCount,sizeX
          call Driver_abortFlash("doInterpExtrap: block index confusion!")
       end if
       if (holeCount .NE. ih) then
          print*,"doInterpExtrap: block index confusion; holeCount,ih=",holeCount,ih
          call Driver_abortFlash("doInterpExtrap: block index confusion!")
       end if
#ifdef DEBUG_MIDVERBOSE
       print*,'Inter/Extra UNK about to fill data holes at',holeInd(1:holeCount),'blkId',blockId
#endif

       iLoc = NGUARD
       do ih = 1,holeCount
          call ut_hunt(xInterp,nInterp,xHole(ih),iLoc)
#ifdef DEBUG_DRIVER
          print*,' ut_hunt for hole #',ih,' returned iLoc=',iLoc
#endif
          if (iLoc .LE. 0 .OR. iLoc .GE. nInterp) then
             call Driver_abortFlash("doInterpExtrap: Invalid iLoc from ut_hunt!")
          end if
          do ivar = 1,NUNK_VARS
             if (ivar  .NE. iKern .AND. &
                  ivar .NE. PDEN_VAR .AND. &
                  ivar .NE. NUP0_VAR .AND. &
                  ivar .NE. NUP1_VAR .AND. &
                  ivar .NE. NUMP_VAR .AND. &
                  .NOT. (sim_unkCellWeight(ivar)>0.0) ) then
                iInterp = 0
                do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
                   if (indInterp(i) == 0) then
                      iInterp = iInterp + 1
                      yInterp(iInterp) = solnData(ivar,i,j,k)
                   end if
                end do
                oldVal = solnData(ivar,holeInd(ih),j,k)
#ifdef DEBUG_DRIVER
99              format('ut_polint for: ',1P,G10.3,'->',G10.3,',',G10.3,'-> ??',G10.3,'??,',G10.3,'->',G10.3,',')
                print 99,xInterp(iLoc),yInterp(iLoc),xHole(ih),oldVal,xInterp(iLoc+1),yInterp(iLoc+1)
#endif
                call ut_polint(xInterp(iLoc),yInterp(iLoc),2,xHole(ih),newVal,dummydy)
#ifdef DEBUG_DRIVER
                print*,'Inter/Extra UNK var #',ivar,'hole#',ih,holeInd(ih),' at',xHole(ih),',old/newVal=',oldVal,newVal,'blkId',blockId
#endif
                solnData(ivar,holeInd(ih),j,k) = newVal

             end if
          end do
          if (ANY(solnData(NUMP_VAR,1:holeInd(ih)-1,j,k).NE.0.0) .AND. &
              ANY(solnData(NUMP_VAR,holeInd(ih)+1:sizeX,j,k).NE.0.0))then
             solnData(NUMP_VAR,holeInd(ih),j,k) = solnData(NUMP_VAR,holeInd(ih),j,k) + 1.0e-30
          end if
       end do
       deallocate(xHole)
       deallocate(holeInd)

    else if (axis==JAXIS) then

       ! !!!!!!!!!!!!!!!!!!!!!!!!!! NOW the Y direction !!!!!!!!! (No Z currently!)

       iInterp = 0
       holeCount = 0
       ! currently first interpolates/extrapolates in X direction...
       do i = blkLimitsGC(LOW,JAXIS)+1, blkLimitsGC(HIGH,JAXIS)-1
          if (solnData(NUMP_VAR,j,i,k) == 0.0) then
             holeCount = holeCount + 1
          end if
       end do
       if (holeCount==0) return
       allocate(holeInd(holeCount))
       allocate(xHole(holeCount))
       holeCount = 0
       do i = blkLimitsGC(LOW,JAXIS)+1, blkLimitsGC(HIGH,JAXIS)-1
          if (solnData(NUMP_VAR,j,i,k) == 0.0) then
             holeCount = holeCount + 1
             holeInd(holeCount) = i
          end if
       end do
       indInterp(:) = 0
       do ih = 1,holeCount
          indInterp(holeInd(ih)) = ih
       end do

       sizeX = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
       call Grid_getCellCoords(JAXIS, blockId, CENTER, .TRUE., xInterp, sizeX)

       iInterp = 0; ih = 0
       do i = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
          if (indInterp(i) == 0) then
             iInterp = iInterp + 1
             if (iInterp< i) xInterp(iInterp) = xInterp(i)
          else  !  (indInterp(i) > 0)
             ih = ih + 1
             xHole(ih) = xInterp(i)
          end if
       end do
       nInterp = iInterp
       if (nInterp + holeCount .NE. sizeX) then
          print*,"doInterpExtrap: block index confusion; nInterp,holeCount,sizeX=",nInterp,holeCount,sizeX
          call Driver_abortFlash("doInterpExtrap: block index confusion!")
       end if
       if (holeCount .NE. ih) then
          print*,"doInterpExtrap: block index confusion; holeCount,ih=",holeCount,ih
          call Driver_abortFlash("doInterpExtrap: block index confusion!")
       end if
#ifdef DEBUG_MIDVERBOSE
       print*,'Inter/Extra UNK about to fill data holes at',holeInd(1:holeCount),'blkId',blockId
#endif

       iLoc = NGUARD
       do ih = 1,holeCount
          call ut_hunt(xInterp,nInterp,xHole(ih),iLoc)
#ifdef DEBUG_DRIVER
          print*,' ut_hunt for hole #',ih,' returned iLoc=',iLoc
#endif
          if (iLoc .LE. 0 .OR. iLoc .GE. nInterp) then
             call Driver_abortFlash("doInterpExtrap: Invalid iLoc from ut_hunt!")
          end if
          do ivar = 1,NUNK_VARS
             if (ivar  .NE. iKern .AND. &
                  ivar .NE. PDEN_VAR .AND. &
                  ivar .NE. NUP0_VAR .AND. &
                  ivar .NE. NUP1_VAR .AND. &
                  ivar .NE. NUMP_VAR .AND. &
                  .NOT. (sim_unkCellWeight(ivar)>0.0) ) then
                iInterp = 0
                do i = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
                   if (indInterp(i) == 0) then
                      iInterp = iInterp + 1
                      yInterp(iInterp) = solnData(ivar,j,i,k)
                   end if
                end do
                oldVal = solnData(ivar,j,holeInd(ih),k)
#ifdef DEBUG_DRIVER
99              format('ut_polint for: ',1P,G10.3,'->',G10.3,',',G10.3,'-> ??',G10.3,'??,',G10.3,'->',G10.3,',')
                print 99,xInterp(iLoc),yInterp(iLoc),xHole(ih),oldVal,xInterp(iLoc+1),yInterp(iLoc+1)
#endif
                call ut_polint(xInterp(iLoc),yInterp(iLoc),2,xHole(ih),newVal,dummydy)
#ifdef DEBUG_DRIVER
                print*,'Inter/Extra UNK var #',ivar,'hole#',ih,holeInd(ih),' at',xHole(ih),',old/newVal=',oldVal,newVal,'blkId',blockId
#endif
                solnData(ivar,j,holeInd(ih),k) = newVal

             end if
          end do
       end do
       if (holeCount==blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)-1 .AND. &
            solnData(NUMP_VAR,j,blkLimitsGC(LOW,JAXIS)+1,k) == 0.0 .AND. &
            solnData(NUMP_VAR,j,blkLimitsGC(HIGH,JAXIS)-1,k) == 0.0 .AND. &
            ALL(solnData(NUMP_VAR,j,blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),k)==0.0) .AND. &
            ALL(solnData(SPECIES_BEGIN:SPECIES_END,j,blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),k)==0.0) .AND. &
            ALL(sim_unkCellWeight(SPECIES_BEGIN:SPECIES_END)==0.0) ) then
          where (solnData(DENS_VAR,j,blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),k) < sim_abundanceFixupMaxDens)
             solnData(C_SPEC,j,blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),k) = 0.5
             solnData(O_SPEC,j,blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),k) = 0.5
          end where
       end if
       deallocate(xHole)
       deallocate(holeInd)
    end if


  end subroutine doInterpExtrap



  subroutine doWeight(blockId)
  
    use Driver_data,         ONLY : dr_globalMe
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get

    integer, intent(in) :: blockId
    integer :: numCellBlk, numCellGC
    integer :: sMethod
    integer, dimension (3) :: bcl, bcr
    integer :: i, j, k, ih, holeCountGC, holeCountBlk, partCountGC, partCountBlk
    integer, allocatable :: holeI(:), holeJ(:), holeK(:)
    integer, allocatable :: holeGCI(:), holeGCJ(:), holeGCK(:)  !gp use
    integer, allocatable :: partGCI(:), partGCJ(:), partGCK(:)  !gp use
    integer :: nummax, wMaxLength, filledCount, fCountSector, n
    integer, allocatable :: filledCellI(:), filledCellJ(:), filledCellK(:)
    integer :: m, ndis, nlev, ic, jc, kc
    real :: r2Cell, mCell, mr2Cell, sumCell 

    integer :: fCountLev
    real, allocatable :: filledCellTh(:)
    real :: thTmp, thMin, thCell
    integer :: indTmp, ntheta
    real :: thetaSearchMin, thetaSearchMax

    real :: massmin,massmax
    integer :: coordSizeX,coordSizeY,coordSizeZ
    real, allocatable :: xCoord(:),yCoord(:),zCoord(:)

    logical :: sector, thSearch, emptyCell=.TRUE., searchCond

    real, allocatable :: soln2(:,:,:,:)

    integer :: ndimen, nparticles, nproperties, ncells

    real(8), allocatable :: particle_locations(:, :)
    real(8), allocatable :: particle_properties(:, :)
    real(8), allocatable :: cell_locations(:, :)
    real(8), allocatable :: mean(:), amplitude(:)
    real(8), allocatable :: interpolant(:, :, :)
    real(8) :: gp_sigma=1.0

    call RuntimeParameters_get("sMethod", sMethod)

    ! backup unweighted data to soln2
    allocate(soln2(NUNK_VARS,blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)))
    do ivar = 1,NUNK_VARS
       do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
       do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
       do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
          soln2(ivar,i,j,k)=solnData(ivar,i,j,k)
       end do
       end do
       end do
    end do

       ! block (+gc) boundary
       bcl(:)=1
       bcr(1)=NXB+2*NGUARD
       bcr(2)=NYB+2*NGUARD
       bcr(3)=NZB+2*NGUARD

       select case (NDIM)
       case(1)
          numCellBlk=NXB
          numCellGC=bcr(1)
       case(2)
          numCellBlk=NXB*NYB
          numCellGC=bcr(1)*bcr(2)
       case(3)
          numCellBlk=NXB*NYB*NZB
          numCellGC=bcr(1)*bcr(2)*bcr(3)
       end select

       ! Search for the number of "holes" in the block of nxb*nyb cells
       holeCountGC = 0
       partCountGC = 0
       holeCountBlk =0
       partCountBlk =0

       do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
       do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
       do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
          if(emptyCell) then
             searchCond=(nint(solnData(NUMP_VAR,i,j,k)) == 0.0)
          else
             searchCond=(nint(solnData(NUMP_VAR,i,j,k))<= 1)
          end if
          if (SearchCond) then
             holeCountGC = holeCountGC + 1
          else
             partCountGC = partCountGC + 1
          end if
       end do
       end do
       end do

       do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
       do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
       do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
          if(emptyCell) then
             searchCond=(nint(solnData(NUMP_VAR,i,j,k)) == 0.0)
          else
             searchCond=(nint(solnData(NUMP_VAR,i,j,k))<= 1)
          end if
          if (SearchCond) then
             holeCountBlk = holeCountBlk + 1
          else
             partCountBlk = partCountBlk + 1
          end if
       end do
       end do
       end do

       ! print block(+GC) information

       if((partCountGC+holeCountGC)/=numCellGC) then 
          write(*,'(A30,3i5)') 'Warning, pGC+hGC/=numCellGC',partCountGC,holeCountGC,numCellGC
       end if

       if((partCountBlk+holeCountBlk)/=numCellBlk) then 
          write(*,'(A30,3i5)') 'Warning, pBlk+hBlk/=numCellBlk',partCountBlk,holeCountBlk,numCellBlk
       end if

       write(*,'(A38,6i5)') 'dr_globalMe,blockId,hBlk,pBlk,hGC,pGC',dr_globalMe,blockId,holeCountBlk,holeCountGC

       ! return if the main block (noGC) is full
       if (holeCountBlk==0) then
          write(*,'(A30,4i5)') 'dr_globalMe,blockId,hBlk=0,hGC',dr_globalMe,blockId,holeCountBlk,holeCountGC
          return
       end if

       ! leave it if the block(+GC) is totally empty, no data to interpolate, return
       if(holeCountBlk==numCellBlk.and.holeCountGC==numCellGC) then
          write(*,'(A38,3i5)') 'totally empty in blockID,hblk,hGC', blockId,holeCountBlk,holeCountGC
          return
       end if

       ! initiate arrays for empty cells (holes) and filled cells (particles)
       if(sim_doWAvg) then
          allocate(holeI(holeCountBlk))
          allocate(holeJ(holeCountBlk))
          allocate(holeK(holeCountBlk))
          holeCountBlk = 0

          do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
          do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
          do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
             if(emptyCell) then
                searchCond=(nint(solnData(NUMP_VAR,i,j,k)) == 0.0)
             else
                searchCond=(nint(solnData(NUMP_VAR,i,j,k))<= 1)   
             end if
             if (searchCond) then
                holeCountBlk = holeCountBlk + 1
                holeI(holeCountBlk) = i
                holeJ(holeCountBlk) = j
                holeK(holeCountBlk) = k
             end if
          end do
          end do
          end do

       elseif(sim_doGP) then
          allocate(holeGCI(holeCountGC))
          allocate(holeGCJ(holeCountGC))
          allocate(holeGCK(holeCountGC))
          holeCountGC = 0

          allocate(partGCI(partCountGC))
          allocate(partGCJ(partCountGC))
          allocate(partGCK(partCountGC))
          partCountGC = 0

          do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
          do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
          do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
             if(emptyCell) then
                searchCond=(nint(solnData(NUMP_VAR,i,j,k)) == 0.0)
             else
                searchCond=(nint(solnData(NUMP_VAR,i,j,k))<= 1)
             end if
             if (searchCond) then
                holeCountGC = holeCountGC + 1
                holeGCI(holeCountGC) = i
                holeGCJ(holeCountGC) = j
                holeGCK(holeCountGC) = k
             else
                partCountGC = partCountGC + 1
                partGCI(partCountGC) = i
                partGCJ(partCountGC) = j
                partGCK(partCountGC) = k
             end if
          end do
          end do
          end do
       else
          print*, 'Warning: an interpolation method must be selected!'
          return
       end if

       ! apply GP

       if(sim_doGP) then
          nparticles=partCountGC
          ncells=holeCountGC
          nproperties=NUNK_VARS
          ndimen=NDIM

          if(nparticles==0) then
             write(*,'(A41,4i5)') 'Warning: dr_globalMe,blockId,nparticles=0,ncells',dr_globalMe,blockId,nparticles,ncells
          end if

          allocate(particle_locations(ndimen, nparticles))
          allocate(particle_properties(nproperties, nparticles))
          allocate(cell_locations(ndimen, ncells))
          allocate(mean(nproperties))
          allocate(amplitude(nproperties))
          allocate(interpolant(2, nproperties, ncells))

          do n=1, nparticles
             i=partGCI(n)
             j=partGCJ(n)
             k=partGCK(n)
             particle_locations(1, n)=i
             particle_locations(2, n)=j
             do ivar = 1,NUNK_VARS
                ! if (ivar .NE. iKern .AND. &
                !    ivar .NE. PDEN_VAR .AND. &
                !    ivar .NE. NUP0_VAR .AND. &
                !    ivar .NE. NUP1_VAR .AND. &
                !    ivar .NE. NUMP_VAR .AND. &
                !    .NOT. (sim_unkCellWeight(ivar)>0.0) ) then
                        particle_properties(ivar, n)=solnData(ivar,i,j,k)
                ! end if
             end do
          end do

          do n=1, ncells
             i=holeGCI(n)
             j=holeGCJ(n)
             k=holeGCK(n)
             cell_locations(1,n)=i
             cell_locations(2,n)=j
          end do

          call gp_interpolation(nparticles, ncells, ndimen, nproperties,  &
                  particle_locations, particle_properties, &
                  cell_locations, gp_sigma, mean, amplitude, interpolant)

          do n=1, ncells
             i=holeGCI(n)
             j=holeGCJ(n)
             k=holeGCK(n)
             do ivar = 1,NUNK_VARS
                if (ivar .NE. iKern .AND. &
                    ivar .NE. PDEN_VAR .AND. &
                    ivar .NE. NUP0_VAR .AND. &
                    ivar .NE. NUP1_VAR .AND. &
                    ivar .NE. NUMP_VAR .AND. &
                    .NOT. (sim_unkCellWeight(ivar)>0.0) ) then
                       !print*, 'after gp_interpolation',interpolant(1,ivar,ncells)
                        soln2(ivar,i,j,k)=interpolant(1, ivar, n)
                end if
             end do  
          end do

          deallocate(particle_locations)
          deallocate(particle_properties)
          deallocate(cell_locations)
          deallocate(mean)
          deallocate(amplitude)
          deallocate(interpolant)

          write(*, '(A50,4i5)') 'end of sim_doGP,blocId,nproperties,nparticles,ncells',blockId,nproperties,nparticles,ncells

       end if ! sim_doGP

       ! apply weighted average
       ! process holes in blk one by one: (1) search; (2) weight

       if(sim_doWAvg) then 
       do ih = 1, holeCountBlk
          i=holeI(ih)
          j=holeJ(ih)
          k=holeK(ih)

          select case (NDIM)
          case(1)
            ! search filled cells in each sector

            nummax=bcr(1)
            allocate(filledCellI(nummax))
            filledCellI(:)=i
            filledCount=0  ! total counted cells for weighting
            wMaxLength=bcr(1)

            do m=1,2
              fCountSector=0  ! 
              searchInlev1D: do nlev=1, wMaxLength
                do ic = i-nlev, i+nlev
                  if(ic>=bcl(1).and.ic<=bcr(1)) then
                     if(abs(ic-i)==nlev) then

                       select case(m)
                       case(1)
                          sector=(ic>i)
                       case(2)
                          sector=(ic<i)
                       end select

                       if(sector) then
                         if((nint(solnData(NUMP_VAR,ic,j,k)) >= 1)) then
                            fCountSector=fCountSector+1
                            filledCount=filledCount+1
                            filledCellI(filledCount)=ic
                         end if ! cells filled
                       end if ! cells in sector

                     end if ! cells in level
                  end if ! cells in area 
                end do
                ! if no pts found, then go to next level
                if(fCountSector>0) exit searchInlev1D
              end do searchInlev1D
            end do !m

            ! do weighting for ALL mass fractions, when we have all points
            do ivar = 1,NUNK_VARS
               if (ivar .NE. iKern .AND. &
                  ivar .NE. PDEN_VAR .AND. &
                  ivar .NE. NUP0_VAR .AND. &
                  ivar .NE. NUP1_VAR .AND. &
                  ivar .NE. NUMP_VAR .AND. &
                  .NOT. (sim_unkCellWeight(ivar)>0.0) ) then

                  mCell = 0
                  r2Cell = 0
                  mr2Cell = 0
                  sumCell = 0
                  ! Weight S_hole by using nearby nlev levels cells: ic
                  ! S_hole = Sum(S_ic * m_ic * r_ic ^ -2) / Sum (m_ic * r_ic ^ -2)

                  do n=1,filledCount
                    ic=filledCellI(n)
                    mCell = SolnData(DENS_VAR,ic,j,k)*((ic+0.5)**3-(ic-0.5)**3)
                    r2Cell = (ic-i)**2
                    sumCell = sumCell+SolnData(ivar,ic,j,k)*mCell/r2Cell
                    mr2Cell = mr2Cell+mCell/r2Cell
                  end do
             
                  if(mr2Cell==0) then
                     soln2(ivar,i,j,k)= 0.0
                  else
                     ! weight cells with 1pt by avg "imaginary pt" with real pt 
                     if(nint(solnData(NUMP_VAR,i,j,k))==0) then 
                           soln2(ivar,i,j,k)= sumCell/mr2Cell
                     else
                           soln2(ivar,i,j,k)=(solnData(ivar,i,j,k)+sumCell/mr2Cell)*0.5
                     end if
                  end if

               end if  ! ivar
            end do   ! ivar  
               
            deallocate(filledCellI) 
  
          case(2)
         
            ! search filled cells in each sector
   
            nummax=bcr(1)*bcr(2)

            allocate(filledCellI(nummax))
            allocate(filledCellJ(nummax))
            allocate(filledCellTh(nummax))

            filledCellI(:)=i
            filledCellJ(:)=j
            filledCellTh(:)=0

            filledCount=0  ! total counted cells for weighting

            wMaxLength=max(bcr(1), bcr(2))

            ! 4 pts searching
            if(sMethod==4) then                  
            do m=1,4
              fCountSector=0  ! 
              searchInlev2D: do nlev=1, wMaxLength
                do ic = i-nlev, i+nlev
                do jc = j-nlev, j+nlev
                  if(ic>=bcl(1).and.ic<=bcr(1).and.jc>=bcl(2).and.jc<=bcr(2)) then
!                     if(abs(ic-i)==nlev.or.abs(jc-j)==nlev) then
                     ndis=(ic-i)**2+(jc-j)**2
                     if(ndis>=nlev**2.and.ndis<(nlev+1)**2) then
 
                       select case(m)
                       case(1)
                          sector=(ic>i.and.jc>=j)
                       case(2)
                          sector=(ic<=i.and.jc>j)
                       case(3)
                          sector=(ic<i.and.jc<=j)
                       case(4)
                          sector=(ic>=i.and.jc<j)
                       end select

                       if(sector) then 
                         if((nint(solnData(NUMP_VAR,ic,jc,k)) >= 1)) then
                            fCountSector=fCountSector+1
                            filledCount=filledCount+1
                            filledCellI(filledCount)=ic
                            filledCellJ(filledCount)=jc
                         end if ! cells filled
                       end if ! cells in sector

                     end if ! cells in level

                  end if ! cells in area 
                end do
                end do
                ! if no pts found, then go to next level
                if(fCountSector>0) exit searchInlev2D
              end do searchInlev2D
            end do !m
            end if ! sMethod

            ! 3pts searching
            if(sMethod==3) then

              thetaSearchMin=0.0
              thetaSearchMax=2*PI

              searchlev: do nlev=1, wMaxLength
                  
              fCountLev=0  ! in one specific searching level 
              do ic = i-nlev, i+nlev
              do jc = j-nlev, j+nlev
                 if(ic>=bcl(1).and.ic<=bcr(1).and.jc>=bcl(2).and.jc<=bcr(2)) then
                    if(abs(ic-i)==nlev.or.abs(jc-j)==nlev) then
                       ! get a list of theta [0, 2PI)
                       thCell=atan2(1.0*(jc-j), 1.0*(ic-i))
                       if(thCell<0) thCell=2*PI+thCell

                       ! skip when thCell is NOT in searching angles

                       if(thetaSearchMin<thetaSearchMax) then
                         thSearch=(thCell>thetaSearchMin.and.thCell<thetaSearchMax)
                       else
                         ! notice .or. it's equivalent to th>Min.and.th<2PI.or.th>=0.and.th<Max
                         thSearch=(thCell>thetaSearchMin.or.thCell<thetaSearchMax)
                       end if

                       if(filledCount<=1) thSearch=.TRUE.

                       if(thSearch) then 
                          if((nint(solnData(NUMP_VAR,ic,jc,k)) >= 1)) then
                              fCountLev=fCountLev+1
                              filledCount=filledCount+1
                              filledCellI(filledCount)=ic
                              filledCellJ(filledCount)=jc
                              filledCellTh(filledCount)=thCell
                          end if ! cells filled
                       endif ! cells are between angles

                    end if  ! cells in level
                 end if ! cells in area
               end do
               end do
            
               ! if no pts found, then go to next level, keep searching angles same
               if(fCountLev==0) cycle searchlev

               ! if total number of pts is 1, then no need do sorting, go to next level
               ! but we need its positions
               if(filledCount<=1) cycle searchlev
                  
               ! sort all chosen pts with selection sort method
               thTmp=0
               indTmp=0
               ! do n=filledCount-fCountLev+1, filledCount
               do n=1, filledCount-1
                  thMin=filledCellTh(n)
                  do m = n+1, filledCount
                     if (filledCellTh(m)<thMin) then
                        thTmp=filledCellTh(m)
                        filledCellTh(m)=filledCellTh(n)
                        filledCellTh(n)=thTmp
                        thMin=filledCellTh(n)

                        indTmp=filledCellI(m)
                        filledCellI(m)=filledCellI(n)
                        filledCellI(n)=indTmp

                        indTmp=filledCellJ(m)
                        filledCellJ(m)=filledCellJ(n)
                        filledCellJ(n)=indTmp
                     end if
                  end do      
               end do

               ntheta=1
               do n=1, filledCount-1
                  if(filledCellTh(n+1)-filledCellTh(n)>1.0e-7)  ntheta=ntheta+1
               end do

               ! determine searching area (angles) for next level

               if(ntheta==1) then
                  cycle searchlev
               elseif(ntheta==2) then
                  if(filledCellTh(filledCount)-filledCellTh(1)-PI<1.0e-7) then
                     thetaSearchMin=0.0
                     thetaSearchMax=2*PI
                  elseif(filledCellTh(filledCount)-filledCellTh(1)>PI) then
                     thetaSearchMin=filledCellTh(1)
                     thetaSearchMax=filledCellTh(filledCount)
                  else
                     thetaSearchMin=filledCellTh(filledCount)
                     thetaSearchMax=filledCellTh(1)
                  end if
                  cycle searchlev
               else
                  do n=1,filledCount-1
                     if(filledCellTh(n+1)-filledCellTh(n)>=PI) then
                        thetaSearchMin=filledCellTh(n)
                        thetaSearchMax=filledCellTh(n+1)
                        cycle searchlev
                     end if
                  end do
                  if(filledCellTh(1)+2*PI-filledCellTh(filledCount)>=PI) then
                     thetaSearchMin=filledCellTh(filledCount)
                     thetaSearchMax=filledCellTh(1)
                  else
                     exit searchlev
                  end if
               end if

            end do searchlev

            end if ! sMethod

            ! do weighting for all mass fractions, when we all weighting points
            do ivar = 1,NUNK_VARS
               if (ivar .NE. iKern .AND. &
                  ivar .NE. PDEN_VAR .AND. &
                  ivar .NE. NUP0_VAR .AND. &
                  ivar .NE. NUP1_VAR .AND. &
                  ivar .NE. NUMP_VAR .AND. &
                  .NOT. (sim_unkCellWeight(ivar)>0.0) ) then

                  mCell = 0
                  r2Cell = 0
                  mr2Cell = 0
                  sumCell = 0
                  ! Weight S_hole by using nearby wlev levels cells: (ic,jc)
                  ! S_hole = Sum(S_ij * m_icjc * r_icjc ^ -2) / Sum (m_icjc * r_icjc ^ -2)

                  do n=1,filledCount
                    ic=filledCellI(n)
                    jc=filledCellJ(n)

                    mCell = SolnData(DENS_VAR,ic,jc,k)*((ic+0.5)**2-(ic-0.5)**2)
                    r2Cell = (ic-i)**2+(jc-j)**2

                    sumCell = sumCell+SolnData(ivar,ic,jc,k)*mCell/r2Cell
                    mr2Cell = mr2Cell+mCell/r2Cell
                  end do

                  if(mr2Cell==0) then
                     soln2(ivar,i,j,k)= 0.0
                  else
                     ! weight cells with 1pt by avg "imaginary pt" with real pt 
                     if(nint(solnData(NUMP_VAR,i,j,k))==0) then 
                           soln2(ivar,i,j,k)= sumCell/mr2Cell
                     else
                           soln2(ivar,i,j,k)=(solnData(ivar,i,j,k)+sumCell/mr2Cell)*0.5
                     end if
                  end if

               end if  ! ivar
            end do   ! ivar


            ! print min/max of mass for surrounding cells
            if(.FALSE.) then
               mCell = 0
               massmax=0.0
               massmin=1.0e40

               do n=1, filledCount
                  ic=filledCellI(n)
                  jc=filledCellJ(n)

                  mCell = SolnData(DENS_VAR,ic,jc,k)*((ic+0.5)**2-(ic-0.5)**2)
                  massmin=min(massmin, mCell)
                  massmax=max(massmax, mCell)
               end do

               coordSizeX = blkLimitsGC(HIGH, IAXIS)
               coordSizeY = blkLimitsGC(HIGH, JAXIS)
               coordSizeZ = blkLimitsGC(HIGH, KAXIS)
               allocate(xCoord(coordSizeX))
               allocate(yCoord(coordSizeY))
               allocate(zCoord(coordSizeZ))
               xCoord=0
               yCoord=0
               zCoord=0
               call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., xCoord, coordSizeX)
               call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., yCoord, coordSizeY)
               call Grid_getCellCoords(KAXIS, blockId, CENTER, .true., zCoord, coordSizeZ)
                
               write(*,'(A10,3e14.6)') 'mass ratio', xCoord(i),yCoord(j),massmin/massmax

               deallocate(xCoord)
               deallocate(yCoord)
               deallocate(zCoord)
            end if

            deallocate(filledCellI)
            deallocate(filledCellJ)
            deallocate(filledCellTh)

          case(3)
             
            ! search first
            nummax=bcr(1)*bcr(2)*bcr(3)

            allocate(filledCellI(nummax))
            allocate(filledCellJ(nummax))
            allocate(filledCellK(nummax))

            filledCellI(:)=i
            filledCellJ(:)=j
            filledCellK(:)=k

            filledCount=0  ! total counted cells for weighting

            wMaxLength=max(bcr(1),bcr(2),bcr(3))

            do m=1,8
              fCountSector=0  ! 
              searchInlev3D: do nlev=1, wMaxLength
                do ic = i-nlev, i+nlev
                do jc = j-nlev, j+nlev
                do kc = k-nlev, k+nlev
                  if(ic>=bcl(1).and.ic<=bcr(1).and.jc>=bcl(2).and.jc<=bcr(2).and.kc>=bcl(3).and.kc<=bcr(3)) then
                     if(abs(ic-i)==nlev.or.abs(jc-i)==nlev.or.abs(kc-k)==nlev) then

                       select case(8)
                       case(1)
                          sector=(ic>i.and.jc>=j.and.kc>=k)
                       case(2)
                          sector=(ic>i.and.jc>=j.and.kc<k)
                       case(3)
                          sector=(ic<=i.and.jc>j.and.kc>=k)
                       case(4)
                          sector=(ic<=i.and.jc>j.and.kc<k)
                       case(5)
                          sector=(ic<i.and.jc<=j.and.kc>=k)
                       case(6)
                          sector=(ic<i.and.jc<=j.and.kc<k)
                       case(7)
                          sector=(ic>=i.and.jc<j.and.kc>=k)
                       case(8)
                          sector=(ic>=i.and.jc<j.and.kc<k)
                       end select

                       if(sector) then
                         if((nint(solnData(NUMP_VAR,ic,jc,k)) >= 1)) then
                            fCountSector=fCountSector+1
                            filledCount=filledCount+1
                            filledCellI(filledCount)=ic
                            filledCellJ(filledCount)=jc
                         end if ! cells filled
                       end if ! cells in sector

                     end if ! cells in level

                  end if ! cells in area 
                end do
                end do
                end do
                ! if no pts found, then go to next level
                if(fCountSector>0) exit searchInlev3D
              end do searchInlev3D
            end do !m
         
             
            ! do weighting for all mass fractions, when we have all points
            do ivar = 1,NUNK_VARS
               if (ivar .NE. iKern .AND. &
                  ivar .NE. PDEN_VAR .AND. &
                  ivar .NE. NUP0_VAR .AND. &
                  ivar .NE. NUP1_VAR .AND. &
                  ivar .NE. NUMP_VAR .AND. &
                  .NOT. (sim_unkCellWeight(ivar)>0.0) ) then

                  mCell = 0
                  r2Cell = 0
                  mr2Cell = 0
                  sumCell = 0
                  ! Weight S_hole by using nearby nlev levels cells: (ic,jc,kc)
                  ! S_hole = Sum(S_icjckc * m_icjckc * r_icjckc ^ -2) / Sum (m_icjckc * r_icjckc ^ -2)
                  do n = 1, filledCount
                     ic=filledCellI(n)
                     jc=filledCellJ(n)
                     kc=filledCellK(n)

                     mCell = SolnData(DENS_VAR,ic,jc,kc) ! if volume elements are identical for all cells
                     r2Cell = (ic-i)**2+(jc-j)**2+(kc-k)**2
                     sumCell = sumCell+SolnData(ivar,ic,jc,kc)*mCell/r2Cell
                     mr2Cell = mr2Cell+mCell/r2Cell
                  end do
             
                  if(mr2Cell==0) then
                     soln2(ivar,i,j,k)= 0.0
                  else
                     ! weight cells with 1pt by avg "imaginary pt" with real pt 
                     if(nint(solnData(NUMP_VAR,i,j,k))==0) then 
                           soln2(ivar,i,j,k)= sumCell/mr2Cell
                     else
                           soln2(ivar,i,j,k)=(solnData(ivar,i,j,k)+sumCell/mr2Cell)*0.5
                     end if
                  end if
                  ! if(mr2Cell==0) then
                  !    SolnData(ivar,i,j,k)= 0.0
                  ! else
                  !    SolnData(ivar,i,j,k)= sumCell/mr2Cell
                  ! end if
                  
               end if  ! ivar
            end do   ! ivar  

            deallocate(filledCellI)
            deallocate(filledCellJ)
            deallocate(filledCellK)
 
          end select

       end do  ! ih=1,holeCountBlk
 
       print*, 'end of sim_doWAvg' 

       end if ! sim_doWAvg

       if(sim_doWAvg) then
          deallocate(holeI)
          deallocate(holeJ)
          deallocate(holeK)
       elseif(sim_doGP) then
          deallocate(holeGCI)
          deallocate(holeGCJ)
          deallocate(holeGCK)
          deallocate(partGCI)
          deallocate(partGCJ)
          deallocate(partGCK)
       end if

       do ivar = 1,NUNK_VARS
          do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
          do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
          do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
             solnData(ivar,i,j,k)=soln2(ivar,i,j,k)
          end do
          end do
          end do
       end do

       deallocate(soln2)
       
  end subroutine doWeight

end subroutine Driver_evolveFlash

