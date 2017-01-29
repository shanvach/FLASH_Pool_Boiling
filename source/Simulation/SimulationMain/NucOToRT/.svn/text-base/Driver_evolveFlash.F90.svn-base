!!****if* source/Simulation/SimulationMain/NucOToRT/Driver_evolveFlash
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
  use Simulation_data, ONLY : sim_doConvolve,sim_doInterpExtrap,sim_doLowerBounds,sim_doEos,&
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

  call Logfile_stamp( 'Entering evolution routine' , '[Driver_evolveFlash]')

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
  

  call Logfile_stamp( 'Entering evolution loop' , '[Driver_evolveFlash]')
  call Timers_start("evolution")


  do dr_nstep = dr_nBegin, dr_nend
     print*, 'do dr_nstep = dr_nBegin, dr_nend ;  ', dr_nstep, dr_nBegin, dr_nend
     
     !!Step forward in time. See bottom of loop for time step calculation.
     call Grid_getLocalNumBlks(localNumBlocks)
     call Grid_getListOfBlocks(LEAF,blockList,blockCount)

     call doLog


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
        else
           print*,'NOT applying interpolation / extrapolation for UNK cells that got no particles'
        end if

!!$     case(4)
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
        else
           print*,'NOT applying lower limits to dens,temp,eint'
        end if

     case(4)



        if (.TRUE.) then
        ! Apply lower bounds on values
           print*,'Computing Sum X_i'
           do lb = 1, blockCount
              call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)
              call Grid_getBlkPtr(blockList(lb), solnData)
              do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                 do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                    do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                       sumZiYi = 0.0
                       solnData(SUMX_MSCALAR,i,j,k) = sum(solnData(SPECIES_BEGIN:SPECIES_END,i,j,k))
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
        

     case(5)
        if (sim_doEos) then
           ! Compute pressure
           print*,'Applying Eos_wrapped, mode',MODE_DENS_TEMP
           do lb = 1, blockCount
              call Eos_wrapped(MODE_DENS_TEMP,blkLimits, blockList(lb))
           end do
        else
           print*,'NOT calling Eos_wrapped'
        end if

     case(6)
        if (ndimOg == 2 .AND. geometryOg == CYLINDRICAL) then
           print*,'Filling VELZ_VAR from VELY_VAR',MODE_DENS_TEMP
           do lb = 1, blockCount
              call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)
              call Grid_getBlkPtr(blockList(lb), solnData)
              do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                 do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                    do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                       solnData(VELZ_VAR,i,j,k) = solnData(VELY_VAR,i,j,k)
                       solnData(VELY_VAR,i,j,k) = 0.0
                    end do
                 end do
              end do
              call Grid_releaseBlkPtr(blockList(lb), solnData)
           end do
        end if
        ! Map to output grid
        call sim_setupOutputGrid(dr_globalMe)
        call sim_mapUnkVarsToOutputGrid
        call sim_shareOutputGrid(dr_globalMe)
     case(7)
        if (dr_globalMe == MASTER_PE) then
           ! write the output grid
           call sim_writeOutputGrid
        end if
     end select


     !output a plotfile before the grid changes
     call Timers_start("IO_output")
     call IO_output( dr_simTime, &
          dr_dt, dr_nstep+1, dr_nbegin, endRun, PLOTFILE_AND_PARTICLEFILE)
     call Timers_stop("IO_output")


!!$     call Timers_start("Grid_updateRefinement")
!!$     call Grid_updateRefinement( dr_nstep, dr_simTime, gridChanged)
!!$     call Timers_stop("Grid_updateRefinement")



     call IO_output(dr_simTime,dr_dt,dr_nstep+1,dr_nbegin,endRun,&
          CHECKPOINT_FILE_ONLY)
     call Timers_stop("IO_output")


  enddo
  !The value of dr_nstep after the loop is (dr_nend + 1) if the loop iterated for
  !the maximum number of times.  However, we need to retain the value that
  !dr_nstep had during the last loop iteration, otherwise the number for nstep
  !that will be stored in a final checkpoint file will be wrong.
  dr_nstep = min(dr_nstep,dr_nend)


  call Timers_stop("evolution")
  call Logfile_stamp( 'Exiting evolution loop' , '[Driver_evolveFlash]')
  if(.NOT.endRun) call IO_outputFinal()

  call Timers_getSummary( max(0,dr_nstep-dr_nbegin+1))


  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")
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

          call Logfile_stamp( strBuff(1:3,:), 3, 2, "step")

       else

          write (numToStr(1:), "(F8.3)") dr_redshift
          write (strBuff(3,1), "(A)") "z"
          write (strBuff(3,2), "(A)") trim(adjustl(NumToStr))

          write (numToStr(1:), "(1PE12.6)") dr_dt
          write (strBuff(4,1), "(A)") "dt"
          write (strBuff(4,2), "(A)") trim(adjustl(NumToStr))

          call Logfile_stamp( strBuff, 4, 2, "step")

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
end subroutine Driver_evolveFlash
