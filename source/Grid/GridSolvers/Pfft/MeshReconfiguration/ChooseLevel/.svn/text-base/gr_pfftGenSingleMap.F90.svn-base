!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/ChooseLevel/gr_pfftGenSingleMap
!!
!! NAME
!!
!!  gr_pfftGenSingleMap
!!
!! SYNOPSIS
!!  
!!  gr_pfftGenSingleMap()
!!
!! DESCRIPTION
!!
!! This routine does not send the actual data, but uses communication
!! to generate the information about the communication pattern to be expected
!! when actual data movement happens.
!!
!!***

subroutine gr_pfftGenSingleMap(solveLevel, leafMapMode)
#include "Flash.h"
#include "constants.h"
#include "Pfft.h"

  use Logfile_interface, ONLY : Logfile_open, Logfile_close, Logfile_stamp
  use gr_pfftData, ONLY : pfft_procGrid, pfft_ndim, pfft_myPE, &
       pfft_me, pfft_commWithTopology, pfft_numProcs, pfft_inLen
  use gr_pfftReconfigData, ONLY : pfft_pencilSize, pfft_procLookup, pfft_maxRefLev
  use gr_pfftinterface, ONLY : gr_pfftGetDestPfftCoords, &
       gr_pfftCreateSendFragment, gr_pfftCommunicateNodeMetaData, gr_pfftGridPointTable
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkCornerID,&
       Grid_getBlkIndexLimits, Grid_getBlkRefineLevel, Grid_getBlkType
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  integer, intent(IN) :: solveLevel
  logical, optional, intent(IN) :: leafMapMode
  integer,dimension(IAXIS:KAXIS) :: actualStride, cornerID, actualBlkSize, &
       blkSize, stride
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer, dimension(MAXBLOCKS) :: blkList
  integer, dimension(1:MDIM) :: numProcs, firstPfftProc, guard, &
       pfftProcCoords, pfftProcOffset

  !Prefix g means global, prefix l means local.
  !gBlkStart: The block start cell in global space at solve level.
  !gBlkEnd: The block end cell in global space at solve level.
  !gBlkFragStart: The block fragment start cell in global space at solve level.
  !gBlkFragEnd: The block fragment end cell in global space at solve level.
  !gPencilStart: The pencil start cell in global space.
  !gPencilEnd: The pencil end cell in global space.
  !lBlkFragStart: The block fragment start cell in local space at solve level.
  !lBlkFragEnd: The block fragment end cell in local space at solve level.
  !lPencilFragStart: The pencil start cell for the FLASH fragment in local space.
  !lPencilFragEnd: The pencil end cell for the FLASH fragment in local space.
  !lBlkFragSize: The block fragment size.
  !lRealBlkFragStart: The block fragment start cell in the "real" block.
  !lRealBlkFragEnd: The block fragment end cell in the "real" block.
  integer, dimension(1:MDIM) :: gBlkStart, gBlkEnd, &
       gBlkFragStart, gBlkFragEnd, gPencilStart, &
       gPencilEnd, lBlkFragStart, lBlkFragEnd, &
       lPencilFragStart, lPencilFragEnd, lBlkFragSize, lActualBlkFragStart, &
       lActualBlkFragEnd, gActualBlkStart

  integer :: axis, blkCount, blockID, fragmentSize, i, iProc, jProc, kProc, &
       totalBlkSize, sumFragmentSize, ierr, pfftProc, dims, procPos, blkType, &
       refLevelFactor, blkRefLev
  logical :: leafMap

  
  !First determine the type of map that we need to generate.
  leafMap = .false.
  if (present(leafMapMode)) then
     if (leafMapMode .eqv. .true.) then
        leafMap = .true.
     end if
  end if

#ifdef FLASH_GRID_PARAMESH
  if (leafMap .eqv. .true.) then
     call Grid_getListOfBlocks(LEAF, blkList, blkCount)
  else
     call Grid_getListOfBlocks(REFINEMENT, blkList, blkCount, solveLevel)
  end if
#else
  call Grid_getListOfBlocks(LEAF, blkList, blkCount) !UG mode.
#endif

  !PFFT grid is specified at "solve" level.  (Note: pfft_maxRefLev=lrefine_max).
  stride(1:MDIM) = 1
  stride(1:NDIM) = 2 ** (pfft_maxRefLev - solveLevel)


  blockLoop: do i = 1, blkCount

     blockID = blkList(i)
     call Grid_getBlkRefineLevel(blockID, blkRefLev)
     call Grid_getBlkType(blockID, blkType)

     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
     actualBlkSize(1:MDIM) = blkLimits(HIGH,1:MDIM) - blkLimits(LOW,1:MDIM) + 1
     guard(1:MDIM) = blkLimits(LOW,1:MDIM) - blkLimitsGC(LOW,1:MDIM)

     call Grid_getBlkCornerID(blockID, cornerID, actualStride)
     
     !Find the size of this block at solve level.
     blkSize(1:MDIM) = 1  !For <3D sims
     blkSize(1:NDIM) = (actualStride(1:NDIM) * actualBlkSize(1:NDIM)) / &
          stride(1:NDIM)
     totalBlkSize = product(blkSize(1:MDIM))

     !Find the start & end cell position of the block at solve level.
     gBlkStart(1:MDIM) = 1  !For <3D sims
     gBlkStart(1:NDIM) = &
          ceiling(real(cornerID(1:NDIM)) / real(stride(1:NDIM)))

     gBlkEnd(1:MDIM) = 1  !For <3D sims
     gBlkEnd(1:NDIM) = gBlkStart(1:NDIM) + blkSize(1:NDIM) - 1



     !We first need to determine how many fragments to create 
     !from this one FLASH block.  We use the numProcs array to 
     !store how many PFFT processors the block is distributed over.
     numProcs(1:MDIM) = 1
     call gr_pfftGetDestPfftCoords(gBlkStart, firstPfftProc)
     
     do axis = 1, pfft_ndim        
        procPos = firstPfftProc(axis)

        fragmentLoop: do 

           if (procPos >= pfft_procGrid(axis)) then 
              call Driver_abortFlash &
                   ("[gr_pfftGenMap]: Going to overrun data structure")
           end if

           gPencilEnd(axis) = & 
                pfft_procLookup(axis) % procInfo(procPos) % globalEndGridPoint

           !If the end position of the FLASH block is beyond the 
           !PFFT processor end position, it means that some of the block must
           !reside on the next PFFT processor.
           if (gPencilEnd(axis) < gBlkEnd(axis)) then          
              procPos = procPos + 1
              numProcs(axis) = numProcs(axis) + 1
           else
              exit
           end if

        end do fragmentLoop
     end do


     !Now we loop over each processor and calculate the fragment size and
     !position.  We then use this information to construct a fragment object 
     !which is then appended to a list.

     !In <3D simulations we must ensure unused dimensions hold the value 1.
     gBlkFragStart = 1; gBlkFragEnd = 1
     lBlkFragStart = 1; lBlkFragEnd = 1
     lPencilFragStart = 1; lPencilFragEnd = 1

     refLevelFactor =  2 ** (abs(blkRefLev - solveLevel)) !+ve integer.
     sumFragmentSize = 0


     kAxisLoop: do kProc = 0, numProcs(KAXIS)-1
        pfftProcCoords(KAXIS) = firstPfftProc(KAXIS) + kProc
        if (kProc == 0) then
           gBlkFragStart(KAXIS) = gBlkStart(KAXIS)
        else
           !There are no "holes" in the pencil grid, so we can calculate 
           !the global fragment start position by adding 1 to the previous 
           !end position.
           gBlkFragStart(KAXIS) = gBlkFragEnd(KAXIS) + 1
        end if


        jAxisLoop: do jProc = 0, numProcs(JAXIS)-1
           pfftProcCoords(JAXIS) = firstPfftProc(JAXIS) + jProc
           if (jProc == 0) then
              gBlkFragStart(JAXIS) = gBlkStart(JAXIS)
           else
              gBlkFragStart(JAXIS) = gBlkFragEnd(JAXIS) + 1
           end if


           iAxisLoop: do iProc = 0, numProcs(IAXIS)-1
              pfftProcCoords(IAXIS) = firstPfftProc(IAXIS) + iProc
              if (iProc == 0) then
                 gBlkFragStart(IAXIS) = gBlkStart(IAXIS)
              else
                 gBlkFragStart(IAXIS) = gBlkFragEnd(IAXIS) + 1
              end if


              do axis = 1, pfft_ndim                  
                   
                 gPencilStart(axis) = pfft_procLookup(axis) &
                      % procInfo(pfftProcCoords(axis)) % globalStartGridPoint
                 gPencilEnd(axis) = pfft_procLookup(axis) &
                      % procInfo(pfftProcCoords(axis)) % globalEndGridPoint                 

                 !The fragment end position is either the pencil end position 
                 !or the block end position.
                 if (gBlkEnd(axis) > gPencilEnd(axis)) then
                    gBlkFragEnd(axis) = gPencilEnd(axis)
                 else
                    gBlkFragEnd(axis) = gBlkEnd(axis)
                 end if

                 !These describe the position of the grid points relative to the
                 !actual FLASH block of interest.
                 lBlkFragStart(axis) = gBlkFragStart(axis) - &
                      gBlkStart(axis) + guard(axis) + 1  
                 lBlkFragEnd(axis) = gBlkFragEnd(axis) - &
                      gBlkStart(axis) + guard(axis) + 1

                 !These describe the position of the grid points relative to the
                 !actual pencil block of interest.
                 lPencilFragStart(axis) = gBlkFragStart(axis) - &
                      gPencilStart(axis) + 1
                 lPencilFragEnd(axis) = lPencilFragStart(axis) + &
                      (gBlkFragEnd(axis) - gBlkFragStart(axis))


                 !Find the local coordinates in the "real" block.
                 !We are operating with respect to the solve refinement level 
                 !and not with respect to the local block's refinement level.
                 !This means we must multiply by refLevelFactor when 
                 !restricting and divide by refLevelFactor when prolonging.
                 if (blkRefLev == solveLevel) then

                    !Same refinement level.  Copy the local coordinates.
                    lActualBlkFragStart(axis) = lBlkFragStart(axis)
                    lActualBlkFragEnd(axis) = lBlkFragEnd(axis)

                 else if (blkRefLev > solveLevel) then

                    !We need to restrict.  Multiply by refLevelFactor.
                    gActualBlkStart = gBlkStart(axis) * refLevelFactor
                    lActualBlkFragStart(axis) = &
                         (gBlkFragStart(axis) * refLevelFactor) - &
                         gActualBlkStart(axis) + guard(axis) + 1
                    lActualBlkFragEnd(axis) = lActualBlkFragStart(axis) + &
                         (lBlkFragSize(axis) * refLevelFactor)

                 else if (blkRefLev < solveLevel) then

                    !We need to prolong.  Divide by refLevelFactor.
                    if (mod(gBlkFragStart(axis),refLevelFactor) /= 0) then
                       !This shouldn't happen unless the user has chosen 
                       !an odd number of grid points for the block size.
                       call Driver_abortFlash &
                            ("Stop before there is an integer rounding error")
                    end if
                    gActualBlkStart(1:MDIM) = gBlkStart(1:MDIM) / refLevelFactor
                    lActualBlkFragStart(axis) = &
                         (gBlkFragStart(axis) / refLevelFactor) - &
                         gActualBlkStart(axis) + guard(axis) + 1
                    lActualBlkFragEnd(axis) = lActualBlkFragStart(axis) + &
                         (lBlkFragSize(axis) / refLevelFactor)
                 end if

              end do

              !Check fragment size is valid.
              fragmentSize = product( & 
                   (gBlkFragEnd(1:NDIM) - gBlkFragStart(1:NDIM) + 1) )
              
              sumFragmentSize = sumFragmentSize + fragmentSize
              if (sumFragmentSize > totalBlkSize) then
                 call Driver_abortFlash ("[gr_pfftGenSingleMap]:" // &
                      "Block fragment calculation error (1)")
              end if

              
              call MPI_Cart_rank(pfft_commWithTopology, pfftProcCoords(1), &
                   pfftProc, ierr)
              
#ifdef DEBUG_PFFT
              print *, "Processor:", pfft_myPE, "sending a fragment to:", &
                   pfftProc, "Start:", gBlkFragStart, &
                   ", end:", gBlkFragEnd
#endif DEBUG_PFFT

              call gr_pfftCreateSendFragment(pfft_myPE, blockID,  &
                   lBlkFragStart, lBlkFragEnd, &
                   lActualBlkFragStart, lActualBlkFragEnd, &
                   blkType, blkRefLev, solveLevel, &
                   pfftProc, lPencilFragStart, lPencilFragEnd)

           end do iAxisLoop
        end do jAxisLoop
     end do kAxisLoop
     !-----------------------------------------------------------------------

     !Check that the sum of fragments sizes is equal to the total block size.
     if (sumFragmentSize /= totalBlkSize) then  
        call Driver_abortFlash("[gr_pfftGenSingleMap]: " // &
             "Block fragment calculation error (2)")
     end if

  end do blockLoop

end subroutine gr_pfftGenSingleMap
