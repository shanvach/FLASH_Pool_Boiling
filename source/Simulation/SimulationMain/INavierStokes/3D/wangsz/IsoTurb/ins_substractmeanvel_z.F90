! subroutine ins_substractmeanvel_z
! Substract mean velocity by slices on Z planes for turbulence simulations:
! --------------------------------------------------------------------------

#include "Flash.h"

  subroutine ins_substractmeanvel_z(OUT_FACE_VAR,blockCount,blockList)

  ! Modules Use:
  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkCenterCoords,&
                             Grid_getBlkPhysicalSize,&
                             Grid_getBlkBoundBox,    &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_getBlkCornerID,    &
                             Grid_getGlobalIndexLimits

  
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY : lrefine_max
  use Grid_data, only : gr_meshMe, gr_meshNumProcs, gr_meshComm, gr_nblockZ
#else
  use Grid_data, only : gr_meshMe, gr_meshNumProcs, gr_meshComm
#endif

  use Driver_data, ONLY : dr_axisNumProcs

  implicit none
#include "constants.h"
#include "Flash_mpi.h"
#include "IncompNS.h"

  integer :: OUT_FACE_VAR
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(blockCount) :: blockList

  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData
  
  real :: del(MDIM),coord(MDIM),bsize(MDIM),boundBox(CONSTANT_TWO,MDIM),dx,dy,dz
  real :: volcell

  ! Total corners:
  integer, save :: TotCorners 
  integer, save, allocatable, dimension(:) :: CorIDList
  real, save, allocatable, dimension(:)    :: volume, volumeaux
  real, save, allocatable, dimension(:,:)  :: velmean, velmeanaux

  integer iID,CorID,k,kk,lb,blockID,ierr

  logical, save :: firstcall = .TRUE.

  integer, dimension(MDIM) :: cornerID,stride

  real, parameter :: eps = 1.e-12 


  ! Number of corners on finer level:
  if (firstcall) then
#ifdef FLASH_GRID_PARAMESH
    TotCorners = gr_nblockZ* 2**(lrefine_max-1) 
#else /* Uniform grid */
    TotCorners = dr_axisNumProcs(KAXIS)
#endif

    ! Allocate mean veloc variables:
    allocate(CorIDList(TotCorners))
    allocate(velmean(TotCorners*NZB,MDIM))
    allocate( volume(TotCorners*NZB))
    allocate(velmeanaux(TotCorners*NZB,MDIM))
    allocate( volumeaux(TotCorners*NZB))


    ! Populate corner ID list for all block sets in Z direction:
    do iID = 1,TotCorners
        CorIDList(iID) = (iID-1)*NZB + 1
    enddo

    firstcall = .false.
    
  endif

  ! Compute meanU, meanV, meanW
  velmean = 0.
  volume  = 0.
  do iID = 1,TotCorners

    CorID = CorIDList(iID)

    do lb = 1,blockCount

        blockID = blockList(lb)
 
        ! Block Corner ID
        call Grid_getBlkCornerID(blockId, cornerID, stride)

        if (cornerID(KAXIS) .ne. CorID) cycle  

        ! Get blocks dx, dy ,dz:
        call Grid_getDeltas(blockID,del)

        ! Cell volume:
        volcell = del(IAXIS)*del(JAXIS)*del(KAXIS) ! dx*dy*dz

        ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
        call Grid_getBlkPtr(blockID,facezData,FACEZ)

        ! Now loop in k:
        kk = CorID - 1
        do k = GRID_KLO , GRID_KHI

          kk = kk + 1

          ! U average:
          velmean(kk,IAXIS) = velmean(kk,IAXIS) + &
          sum(facexData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,k))*volcell

          ! V average:
          velmean(kk,JAXIS) = velmean(kk,JAXIS) + &
          sum(faceyData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,k))*volcell

          ! W average:
          velmean(kk,KAXIS) = velmean(kk,KAXIS) + &
          sum(facezData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,k))*volcell

          ! Sum to processors total volume:
          volume(kk) = volume(kk) + real(NXB*NYB*1)*volcell

        enddo

        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
        
    enddo

  enddo

  ! Now all reduce:
  velmeanaux = velmean
  call MPI_Allreduce(velmeanaux, velmean, (TotCorners*NZB)*MDIM, FLASH_REAL,&
                     MPI_SUM, gr_meshComm, ierr)
  volumeaux  = volume
  call MPI_Allreduce(volumeaux, volume, (TotCorners*NZB), FLASH_REAL,&
                     MPI_SUM, gr_meshComm, ierr)  
  
  ! Do plane volume averages:
  do iID = 1,TotCorners
      do kk= (iID-1)*NZB+1 , iID*NZB

         if (volume(kk) .lt. eps) then
           write(*,*) ' '
           write(*,*) 'Mype=',gr_meshMe,', slice kk=',kk
           write(*,*) 'Plane Volume less than eps, vol=',volume(kk)
           call Driver_abortFlash("Plane volume zero, this abort is to check single lref runs.") 
         endif

         velmean(kk,IAXIS:KAXIS) = velmean(kk,IAXIS:KAXIS)/volume(kk)

      enddo
  enddo  


  ! Finally substract and assign:
  do iID = 1,TotCorners

    CorID = CorIDList(iID)

    do lb = 1,blockCount

        blockID = blockList(lb)
 
        ! Block Corner ID
        call Grid_getBlkCornerID(blockId, cornerID, stride)

        if (cornerID(KAXIS) .ne. CorID) cycle

        ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
        call Grid_getBlkPtr(blockID,facezData,FACEZ)
 
        kk = CorID - 1
        do k = GRID_KLO , GRID_KHI
          kk = kk + 1

          ! U vels:
          facexData(OUT_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,k) = &
          facexData(VELC_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,k) - velmean(kk,IAXIS)         

          ! V vels:
          faceyData(OUT_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,k) = &
          faceyData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,k) - velmean(kk,JAXIS)
          
          ! W vels:
          facezData(OUT_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,k) = &
          facezData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,k) - velmean(kk,KAXIS)
 
 
        enddo

        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

    enddo

  enddo

  return
  end subroutine



