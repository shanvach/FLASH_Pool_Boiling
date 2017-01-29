!!****if* source/physics/SolidMechanics/SolidMechanicsMain/sm_gatherMassCent.F90
!!
!! NAME
!!
!!
!!
!! SYNOPSIS
!!
!!
!! VARIABLES
!!
!!
!! DESCRIPTION
!!
!!  Shizhao Wang modified from the file SolidMechanics_init.F90
!!  Gather the mass centers of the solid bodies to all processors
!!  The mass centers will be used in the collision model.
!!
!!  Nov 05, 2014
!!
!!***

subroutine sm_gatherMassCent()
  use SolidMechanics_data, only : sm_MeshMe,sm_NumProcs, sm_NumBodies, sm_BodyInfo, sm_meshComm, &
                                  sm_metaBodies, sm_pos
  use Driver_interface, only : Driver_getMype,Driver_getNumProcs, &
                               Driver_abortFlash, &
                               Driver_getComm

  implicit none
#include "constants.h"
#include "SolidMechanics.h"
#include "sm_integrator.h"
#include "Flash_mpi.h"
#include "Flash.h"

  integer :: smBodyCount, maxBodyCount
  integer, allocatable, dimension(:) :: smBodyList, bodyCounts
  integer, allocatable, dimension(:,:) :: bodyLists
  real, allocatable, dimension(:,:) :: smMassCenter
  real, allocatable, dimension(:,:,:) :: massCenters

  integer :: iCnt, iProcs, ind, ibd, imaster, ierr, cnt

  if(sm_metaBodies > 0) then
      call Driver_abortFlash('The routine sm_gatherMassCent has not been tested for meta-bodies.')
  endif

  ! Which Processor am I:
  call Driver_getMype(GLOBAL_COMM, sm_MeshMe)

  ! Number of Processors:
  call Driver_getNumProcs(GLOBAL_COMM, sm_NumProcs)

  ! Get Comm
  call Driver_getComm( GLOBAL_COMM, sm_meshComm )

  smBodyCount = 0
  Do ibd = 1, sm_NumBodies
    if(sm_BodyInfo(ibd)%BodyMaster == sm_MeshMe) then
       smBodyCount = smBodyCount + 1
    endif
  Enddo

  allocate(bodyCounts(sm_NumProcs))
  bodyCounts = 0
  call MPI_Allgather(smBodyCount,1,FLASH_INTEGER,bodyCounts(1),1,FLASH_INTEGER,sm_meshComm, ierr)
  
  maxBodyCount = maxval(bodyCounts)

  !debug
!  write(*,*) 'smBodyCount:', smBodyCount, 'id:', sm_meshMe
!  write(*,*) 'BodyCounts:', bodyCounts, 'id:', sm_meshMe
!  write(*,*) 'maxBodyCount:', maxBodyCount, 'id:', sm_meshMe

  allocate(smBodyList(maxBodyCount))
  smBodyList = 0

  allocate(smMassCenter(NDIM,maxBodyCount))
  smMassCenter = 0.
  
  imaster = BODYFRAME_NODE

  smBodyCount = 0
  Do ibd = 1, sm_NumBodies
    if(sm_BodyInfo(ibd)%BodyMaster == sm_MeshMe) then
       smBodyCount = smBodyCount + 1
       smBodyList(smBodyCount) = ibd
       smMassCenter(1,smBodyCount) = sm_BodyInfo(ibd)%x(imaster) & 
                               &   + sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(sm_BodyInfo(ibd)%ix,imaster))
       smMassCenter(2,smBodyCount) = sm_BodyInfo(ibd)%y(imaster) &
                               &   + sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(sm_BodyInfo(ibd)%ix+1,imaster))
#if NDIM == MDIM
       smMassCenter(3,smBodyCount) = sm_BodyInfo(ibd)%z(imaster) &
                               &   + sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(sm_BodyInfo(ibd)%ix+2,imaster))
#endif
    endif
  Enddo

  allocate(bodyLists(maxBodyCount,sm_NumProcs))
  bodyLists = 0

  allocate(massCenters(NDIM,maxBodyCount,sm_NumProcs))
  massCenters = 0

  call MPI_Allgather(smBodyList(1),maxBodyCount,FLASH_INTEGER,bodyLists(1,1),maxBodyCount,FLASH_INTEGER,sm_meshComm, ierr)
 
  cnt = NDIM*maxBodyCount 
  call MPI_Allgather(smMassCenter(1,1),cnt,FLASH_REAL,massCenters(1,1,1),cnt,FLASH_REAL,sm_meshComm, ierr)

  do iProcs = 1, sm_NumProcs
    if(bodyCounts(iProcs)>0) then
      do iCnt = 1, bodyCounts(iProcs)
        ind = bodyLists(iCnt, iProcs)
        sm_pos(:,ind) = massCenters(:,iCnt,iProcs)
      enddo
    endif
  enddo
  
!  write(*,*) 'smbodyLists:', smBodyList, 'id:', sm_meshMe
!  write(*,*) 'smMassCenter:', smMassCenter, 'id:', sm_meshMe
!  write(*,*) 'BodyLists:', bodyLists, 'id:', sm_meshMe
!  write(*,*) 'massCenters:', masscenters, 'id:', sm_meshMe
!  write(*,*) 'sm_pos:', sm_pos, 'id:', sm_meshMe

  deallocate(smBodyList, smMassCenter)
  deallocate(bodyLists, bodyCounts, massCenters)

  return

end subroutine sm_gatherMassCent
