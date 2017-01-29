

! Fills neighProcsCount and neighProcsList of superlu_common

#include "Flash.h"
#include "constants.h"
#include "Superlu.h"

subroutine gr_slugetSubsetNeighProcs(blockList_set,blockCount_set,diagblk_flag)

  use superlu_common, only : neighProcsCount,neighProcsList

  use Grid_data, only : gr_meshMe, gr_meshNumProcs

  use tree, only : surr_blks

  use ut_qsortInterface, ONLY : ut_qsort
  use gr_interfaceTypeDecl

  implicit none
  integer, intent(in) :: blockList_set(MAXBLOCKS),blockCount_set
  integer, intent(in) :: diagblk_flag

  ! Local variables
  integer, parameter :: CONSTANT_THREE = 3
  integer :: ENDI, ENDJ, ENDK, allcenters, procSize, procCount, procID
  integer, allocatable, dimension(:) :: procs
  integer :: i,j,k, b, p, u, n, blockID, numNegh  
  integer :: vec_ij(CONSTANT_TWO),vec_jk(CONSTANT_TWO),vec_ik(CONSTANT_TWO)

  type (AllBlockRegions_t) :: surrBlksSummary
  integer,dimension(BLKNO:TYPENO) :: negh_prop

  logical :: blktest_flg


  ! Find which are the processors that share block boundaries with me:
  ! gives me neighProcsCount  
  ! call gr_pdsluFindNeighs(blockCount_set,blockList_set,...)
  ENDI = CONSTANT_THREE
  ENDJ = max(1,K2D*CONSTANT_THREE)
  ENDK = max(1,K3D*CONSTANT_THREE)
  allCenters = 2**NDIM

  ! Up to "2**(NDIM-1)" neighbors for each guard cell region
  ! Exactly "((ENDI*ENDJ*ENDK)-1)" guard cell regions in each block
  ! Exactly "blockCount_set" blocks in this MPI rank.
  ! Add "1" to ensure my MPI rank is in the procs array.
  procSize = (2**(NDIM-1) * ((ENDI*ENDJ*ENDK)-1) * blockCount_set) + 1
  allocate(procs(procSize))

  procCount = 0
  do b = 1, blockCount_set
     blockID = blockList_set(b)

     ! The call to this routine assumes the block is a leaf block,
     ! It should work for lrefinements that cover the whole computaitonal domain
     ! but will not work for parents.
     do k = 1, ENDK
        do j = 1, ENDJ
           do i = 1, ENDI

           !Initialize all neighbor details in region (i,j,k) to NONEXISTENT.
           surrBlksSummary % regionInfo(i,j,k) % numNegh = NONEXISTENT
           surrBlksSummary % regionInfo(i,j,k) &
                % details(BLKNO:TYPENO, 1:2**(NDIM-1)) = NONEXISTENT

           select case(diagblk_flag)
           case(DIAGBLKS_YES)
              blktest_flg = ((i*j*k) /= allCenters)
           case(DIAGBLKS_NO)
              vec_ij = (/ i , j /)
              vec_jk = (/ j , k /)
              vec_ik = (/ i , k /)
              blktest_flg = ( ((i*j*k) /= allCenters) .and.  &
                                  (all(vec_ij .eq. 2)  .or.  &
                                   all(vec_jk .eq. 2)  .or.  &
                                   all(vec_ik .eq. 2)) ) ! No Edge or corner Blocks needed 1 level compute.
           case default
              call Driver_abortFlash("gr_slugetSubsetNeighProcs : incorrect diagblk_flag.")
           end select

           if ( blktest_flg ) then 

              !May not be a block at this position.  We may be at 
              !the edge of the domain, and so surr_blks at i,j,k contains 
              !the external boundary conditions.
              negh_prop(:) = surr_blks(:,i,j,k,blockID)

              !First check for an external boundary.
              if (negh_prop(BLKNO) <= PARAMESH_PHYSICAL_BOUNDARY) then

                 !This is the case for e.g. an external reflecting boundary.
                 !...Just copy the external boundary conditions.
                 surrBlksSummary % regionInfo(i,j,k) % numNegh = 0
                 surrBlksSummary % regionInfo(i,j,k) &
                      % details(BLKNO:TYPENO,1) = negh_prop(BLKNO:TYPENO)

              !Now check for wacky values in surr_blks.
              else if ( &
                   (negh_prop(BLKNO) < -1) .or. &
                   (negh_prop(BLKNO) == 0) .or. &
                   (negh_prop(BLKNO) > MAXBLOCKS) .or. &
                   (negh_prop(PROCNO) < -1) .or. &
                   (negh_prop(PROCNO) >= gr_meshNumProcs) &
                   ) then

                 print *, gr_meshMe, negh_prop(:)
                 call Driver_abortFlash("Unexpected surr_blks values")


              else

                 !! This is the situation when the neighbor is at the
                 !! same level of resulution. There is only one neighbor
                 !! very simply found.
                 surrBlksSummary % regionInfo(i,j,k) % numNegh = 1
                 surrBlksSummary % regionInfo(i,j,k) &
                      % details(BLKNO:TYPENO,1) = negh_prop(BLKNO:TYPENO)

              endif


              numNegh = surrBlksSummary % regionInfo(i,j,k) % numNegh
              do n = 1, numNegh
                 procID = surrBlksSummary % regionInfo(i,j,k) % &
                          details(PROCNO,n)
                 if (procID /= gr_meshMe) then
                    procCount = procCount + 1
                    procs(procCount) = procID
                 end if
              end do
           end if
           end do
        end do
     end do
  end do

if (allocated(neighProcsList)) deallocate(neighProcsList)
  allocate(neighProcsList(procCount+1))
  if(procCount > 0) then
     call ut_qsort(procs, procCount)
     neighProcsList = -1
     neighProcsList(1) = gr_meshMe
     u = 1
     do p = 1, procCount
        if (.not.(any(procs(p)==neighProcsList))) then
           u = u + 1
           neighProcsList(u) = procs(p)
        end if
     end do
  else
     u = 1
     neighProcsList(u) = gr_meshMe
  end if
  deallocate(procs)
  neighProcsCount = u

  return


end subroutine gr_slugetSubsetNeighProcs
