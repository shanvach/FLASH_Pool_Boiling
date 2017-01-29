
! Creates Chunk of Matrix belonging to Processor in CSR format

#include "Flash.h"
#include "constants.h"

subroutine gr_sluBuildAdist(A, m, n)

  use Grid_data, only : gr_meshComm,gr_meshMe

  use superlu_common, only :  neighProcsCount,neighProcsList,     &
                              neighProc_blkcnt,neighProc_blkList, &
                              neighProc_blkList_idx,              &
                              blockCount_start_idx,               &
                              nloc_dofs,n_glob,m_glob,            &
                              nnz_glob,m_loc,nnz_loc,fst_row,     &
                              nnz_allocate,nzval,colind,rowptr,dofs_block

  use superlu_mod


  use Grid_interface,   ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,   &
                               Grid_getBlkIndexLimits,Grid_getBlkBC, &
                               Grid_getDeltas


  use tree, only : surr_blks

  implicit none
#include "Flash_mpi.h"
  integer(superlu_ptr), intent(in) :: A
  integer, intent(inout) :: m, n
  !integer, intent(inout) :: colind(nnz_allocate), rowptr(m_loc+1)  
  !real, intent(inout) :: nzval(nnz_allocate)

  ! Local Variables:
  integer :: blockList_set(MAXBLOCKS),blockCount_set
  integer,dimension(BLKNO:TYPENO) :: negh_prop
  integer, dimension(2,MDIM):: faces 
  real :: row_fact(2*NDIM+1)

  real :: ax,bx,cx,ay,by,cy,az,bz,cz,nzvali(10),nzvali2(10),bbb
  integer :: nnzloci,nnzglobi,nnz_loci,colindi(10),colindi2(10),blockID,p,i,j,k,ii,jj,kk,iii,jjj
  integer :: idxb, ierr, lb, negh_blk, negh_proc, nghproc, nrowi, aaa
  integer :: idof_im1,idof_ip1,idof_i
  integer :: idof_jm1,idof_jp1,idof_j
  integer :: idof_km1,idof_kp1,idof_k
  real :: ss,uu,pp,ee,rr,ll

  integer,dimension(2,MDIM)::blkLimits,blkLimitsGC
  real, dimension(MDIM)     :: del
  real, pointer, dimension(:,:,:,:) :: facevarx,facevary,facevarz

  character(len=4) :: ind_me

  ! Local number of A rows
  m_loc        =  nloc_dofs
  nnz_allocate = (2*NDIM + 1)*nloc_dofs

  ! Global Numeration for first row in the processor
  fst_row    = dofs_block*blockCount_start_idx     !Zero-based

  ! allocate 
  if(allocated(nzval))  deallocate(nzval)
  if(allocated(colind)) deallocate(colind)
  if(allocated(rowptr)) deallocate(rowptr)
  allocate(nzval(nnz_allocate),colind(nnz_allocate),rowptr(m_loc+1))
     
  ! Mype Blocks
  blockCount_set = neighProc_blkcnt(CONSTANT_TWO,CONSTANT_ONE)      ! subset blocks number in mype
  blockList_set(1:blockCount_set)=neighProc_blkList(1:blockCount_set,CONSTANT_ONE) ! subset blocks in mype 
!  print*,"Proc",gr_meshMe,"FirstRow",fst_row,"Subset",blockCount_set,"-->",blockList_set(1:blockCount_set)

  nnzglobi = 0
  nnz_loc  = 0
  nrowi    = 0
  do lb=1,blockCount_set

     blockID = blockList_set(lb) ! This is in Processors

     ! Get Block BCs:
     call Grid_getBlkBC (blockID, faces)     

     ! Get deltas:
     call Grid_getDeltas(blockID, del)

     ! Get Block limits
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Get Densities:
     call Grid_getBlkPtr(blockID,facevarx ,FACEX)
     call Grid_getBlkPtr(blockID,facevary ,FACEY)
#if NDIM == 3
     call Grid_getBlkPtr(blockID,facevarz ,FACEZ)
#endif     
 
     
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

              ! Set the Row Ptr
              nrowi = nrowi + 1
              rowptr(nrowi) = nnz_loc 

              ! Numeration from 1:NXB, etc.
              ii = i - NGUARD
              jj = j - NGUARD
#if NDIM == MDIM
              kk = k - NGUARD
#else
              kk = 1
#endif

              ! Initialize local row arrays
              nnzloci = 1
              nzvali  = 0.
              colindi = 0

              ! Discretization Coefficients in x
              !ax =  1./(del(IAXIS)**2.0)
              !bx = -2./(del(IAXIS)**2.0) 
              !cx =  1./(del(IAXIS)**2.0)
              ax =  facevarx(MGW8_FACE_VAR,i  ,j  ,k)/(del(IAXIS)**2.0)
              bx = -(facevarx(MGW8_FACE_VAR,i  ,j  ,k)+facevarx(MGW8_FACE_VAR,i+1,j  ,k))/(del(IAXIS)**2.0)
              cx =  facevarx(MGW8_FACE_VAR,i+1,j  ,k)/(del(IAXIS)**2.0)

!print*,"Xden",i,j,k,facevarx(MGW8_FACE_VAR,i  ,j  ,k),facevarx(MGW8_FACE_VAR,i+1,j  ,k)

              ! ax coefficient:
              if( i .eq. blkLimits(LOW,IAXIS)) then

                 select case(faces(LOW,IAXIS))

                 !KPD - 0:3,5 & 1:2,3 - Interior Low X Bdrys
                 case(NOT_BOUNDARY,PERIODIC)

                    ! Find the dof of neighbor block in i-1,j,k
                    negh_proc = surr_blks(PROCNO,LEFT_EDGE,CENTER,CENTER,blockID)                    
                    negh_blk  = surr_blks(BLKNO ,LEFT_EDGE,CENTER,CENTER,blockID)                    

!print*,"Proc",gr_meshMe,lb,blockID,"PERIODIC",negh_proc,negh_blk

                    nghproc = -1
                    do p=1,neighProcsCount
                       if (neighProcsList(p) .eq. negh_proc) nghproc = p
                    enddo

                    ! Dof in global numbering:
                    !KPD - idxb = global block index (within subset) of neighbor
                    !KPD - idof_im1 = global cell index of neighbor
                    idxb = neighProc_blkList_idx(negh_blk,nghproc)
                    idof_im1 = NXB + (jj-1)*NXB + (kk-1)*NXB*NYB + (idxb-1)*dofs_block
!if (gr_meshMe .eq. 0 .AND. blockID .eq. 3) then
!print*,"Proc",gr_meshMe,"blk",blockID,"ijk",i,j,k,"idxb",idxb,"idof_im1",idof_im1
!end if
                    
                    ! Assign to Matrix:
                    nnzglobi         = nnzglobi + 1
                    nnzloci          = nnzloci  + 1    
                    nzvali(nnzloci)  = ax                    
                    colindi(nnzloci) = idof_im1 - 1 !zero-based

!                    !KPD - Incompressible Reference Pressure.
!                    if (lb.eq.1 .AND. i .eq. blkLimits(LOW,IAXIS) .AND. &
!                                      j .eq. blkLimits(LOW,JAXIS) .AND. &
!                                      k .eq. blkLimits(LOW,KAXIS) ) then
!                       nzvali(nnzloci)  = 0.0
!                    end if 
 
                 case(OUTFLOW)

                    ! Add to bx:
                    bx = bx + ax

                 case(DIRICHLET)

                    ! Add to bx:
                    bx = bx - ax   
                 

                 end select


              else

                    ! Dof in global numbering:
                    !KPD - nghproc = 1 ia the index for itself
                    !KPD - Set global block index to itself
                    !KPD - idof_im1 = global cell index of neighbor to left
                    nghproc = 1
                    idxb = neighProc_blkList_idx(blockID,nghproc)
                    idof_im1 = (ii-1) + (jj-1)*NXB + (kk-1)*NXB*NYB + (idxb-1)*dofs_block
                    
!if (gr_meshMe .eq. 0 .AND. blockID .eq. 1) then
!print*,"Proc",gr_meshMe,"blk",blockID,"ijk",i,j,k,"idxb",idxb,"idof_im1",idof_im1
!end if

                    ! Assign to Matrix:
                    nnzglobi         = nnzglobi + 1
                    nnzloci          = nnzloci  + 1   
                    nzvali(nnzloci)  = ax                    
                    colindi(nnzloci) = idof_im1 - 1 !zero-based

!                    !KPD - Incompressible Reference Pressure.
!                    if (lb.eq.1 .AND. i .eq. blkLimits(LOW,IAXIS) .AND. &
!                                      j .eq. blkLimits(LOW,JAXIS) .AND. &
!                                      k .eq. blkLimits(LOW,KAXIS) ) then 
!                       nzvali(nnzloci)  = 0.0
!                    end if


              endif


              ! cx coefficient:
              if ( i .eq. blkLimits(HIGH,IAXIS)) then

                 select case(faces(HIGH,IAXIS))

                 case(NOT_BOUNDARY,PERIODIC)

                    ! Find the dof of neighbor block in i+1,j,k
                    negh_proc = surr_blks(PROCNO,RIGHT_EDGE,CENTER,CENTER,blockID)                    
                    negh_blk  = surr_blks(BLKNO ,RIGHT_EDGE,CENTER,CENTER,blockID)                    

!print*,"Proc",gr_meshMe,lb,blockID,"PERIODIC",negh_proc,negh_blk

                    nghproc = -1
                    do p=1,neighProcsCount
                       if (neighProcsList(p) .eq. negh_proc) nghproc = p
                    enddo

                    ! Dof in global numbering:
                    idxb = neighProc_blkList_idx(negh_blk,nghproc)
                    idof_ip1 = 1 + (jj-1)*NXB + (kk-1)*NXB*NYB + (idxb-1)*dofs_block
                    
!if (gr_meshMe .eq. 0 .AND. blockID .eq. 2) then
!print*,"Proc",gr_meshMe,"blk",blockID,"ijk",i,j,k,"idxb",idxb,"idof_ip1",idof_ip1
!end if


                    ! Assign to Matrix:
                    nnzglobi         = nnzglobi + 1
                    nnzloci          = nnzloci  + 1
                    nzvali(nnzloci)  = cx                    
                    colindi(nnzloci) = idof_ip1 - 1 !zero-based

!                    !KPD - Incompressible Reference Pressure.
!                    if (lb.eq.1 .AND. i .eq. blkLimits(LOW,IAXIS) .AND. &
!                                      j .eq. blkLimits(LOW,JAXIS) .AND. &
!                                      k .eq. blkLimits(LOW,KAXIS) ) then 
!                       nzvali(nnzloci)  = 0.0
!                    end if


                 case(OUTFLOW)

                    ! Add to bx:
                    bx = bx + cx       

                 case(DIRICHLET)

                    ! Add to bx:
                    bx = bx - cx   
              
                 end select

              else

                    ! Dof in global numbering:
                    nghproc = 1
                    idxb = neighProc_blkList_idx(blockID,nghproc)
                    idof_ip1 = (ii+1) + (jj-1)*NXB + (kk-1)*NXB*NYB + (idxb-1)*dofs_block
                    
!if (gr_meshMe .eq. 0 .AND. blockID .eq. 1) then
!print*,"Proc",gr_meshMe,"blk",blockID,"ijk",i,j,k,"idxb",idxb,"idof_ip1",idof_ip1
!end if

                    ! Assign to Matrix:
                    nnzglobi         = nnzglobi + 1
                    nnzloci          = nnzloci  + 1
                    nzvali(nnzloci)  = cx                    
                    colindi(nnzloci) = idof_ip1 - 1 !zero-based

!                    !KPD - Incompressible Reference Pressure.
!                    if (lb.eq.1 .AND. i .eq. blkLimits(LOW,IAXIS) .AND. &
!                                      j .eq. blkLimits(LOW,JAXIS) .AND. &
!                                      k .eq. blkLimits(LOW,KAXIS) ) then 
!                       nzvali(nnzloci)  = 0.0
!                    end if

              end if

              nnzglobi         = nnzglobi + 1

              ! bx coefficient
              nghproc = 1
              idxb = neighProc_blkList_idx(blockID,nghproc)
              idof_i = ii + (jj-1)*NXB + (kk-1)*NXB*NYB + (idxb-1)*dofs_block
                    
!if (gr_meshMe .eq. 0 .AND. blockID .eq. 2) then
!print*,"Proc",gr_meshMe,"blk",blockID,"ijk",i,j,k,"idxb",idxb,"idof_i",idof_i
!end if

              ! Assign to Matrix:
              nzvali(CONSTANT_ONE)  = nzvali(CONSTANT_ONE) + bx                    
              colindi(CONSTANT_ONE) = idof_i - 1 !zero-based              



              ! Discretization coefficients in y:
              !ay = 1./(del(JAXIS)**2.0)
              !by =-2./(del(JAXIS)**2.0)
              !cy = 1./(del(JAXIS)**2.0) 
              ay = facevary(MGW8_FACE_VAR,i  ,j  ,k)/(del(JAXIS)**2.0)
              by =-(facevary(MGW8_FACE_VAR,i  ,j  ,k)+facevary(MGW8_FACE_VAR,i  ,j+1,k))/(del(JAXIS)**2.0)
              cy = facevary(MGW8_FACE_VAR,i  ,j+1,k)/(del(JAXIS)**2.0)



              ! ay coefficient:
              if( j .eq. blkLimits(LOW,JAXIS)) then

                 select case(faces(LOW,JAXIS))

                 case(NOT_BOUNDARY,PERIODIC)

                    ! Find the dof of neighbor block in i,j-1,k
                    negh_proc = surr_blks(PROCNO,CENTER,LEFT_EDGE,CENTER,blockID)                    
                    negh_blk  = surr_blks(BLKNO ,CENTER,LEFT_EDGE,CENTER,blockID)                    

                    nghproc = -1
                    do p=1,neighProcsCount
                       if (neighProcsList(p) .eq. negh_proc) nghproc = p
                    enddo

                    ! Dof in global numbering:
                    idxb = neighProc_blkList_idx(negh_blk,nghproc)
                    idof_jm1 = ii + (NYB-1)*NXB + (kk-1)*NXB*NYB + (idxb-1)*dofs_block

!if (gr_meshMe .eq. 1 ) then!.AND. blockID .eq. 1) then
!print*,"Proc",gr_meshMe,"blk",blockID,"ijk",i,j,k,"idxb",idxb,"idof_ip1",idof_jm1
!end if
                    
                    ! Assign to Matrix:
                    nnzglobi         = nnzglobi + 1
                    nnzloci          = nnzloci  + 1
                    nzvali(nnzloci)  = ay                    
                    colindi(nnzloci) = idof_jm1 - 1 !zero-based

!                    !KPD - Incompressible Reference Pressure.
!                    if (lb.eq.1 .AND. i .eq. blkLimits(LOW,IAXIS) .AND. &
!                                      j .eq. blkLimits(LOW,JAXIS) .AND. &
!                                      k .eq. blkLimits(LOW,KAXIS) ) then
!                       nzvali(nnzloci)  = 0.0
!                    end if
 
                 case(OUTFLOW)

                    ! Add to by:
                    by = by + ay   
 
                 case(DIRICHLET)

                    ! Add to by:
                    by = by - ay                    

                 end select


              else

                    ! Dof in global numbering:
                    nghproc = 1
                    idxb = neighProc_blkList_idx(blockID,nghproc)
                    idof_jm1 = ii + (jj-2)*NXB + (kk-1)*NXB*NYB + (idxb-1)*dofs_block
                    
!if (gr_meshMe .eq. 1 ) then!.AND. blockID .eq. 1) then
!print*,"Proc",gr_meshMe,"blk",blockID,"ijk",i,j,k,"idxb",idxb,"idof_ip1",idof_jm1
!end if
                    
                    ! Assign to Matrix:
                    nnzglobi         = nnzglobi + 1
                    nnzloci          = nnzloci  + 1
                    nzvali(nnzloci)  = ay                    
                    colindi(nnzloci) = idof_jm1 - 1 !zero-based

!                    !KPD - Incompressible Reference Pressure.
!                    if (lb.eq.1 .AND. i .eq. blkLimits(LOW,IAXIS) .AND. &
!                                      j .eq. blkLimits(LOW,JAXIS) .AND. &
!                                      k .eq. blkLimits(LOW,KAXIS) ) then
!                       nzvali(nnzloci)  = 0.0
!                    end if

              endif


              ! cy coefficient:
              if ( j .eq. blkLimits(HIGH,JAXIS)) then

                 select case(faces(HIGH,JAXIS))

                 case(NOT_BOUNDARY,PERIODIC)

                    ! Find the dof of neighbor block in i,j+1,k
                    negh_proc = surr_blks(PROCNO,CENTER,RIGHT_EDGE,CENTER,blockID)                    
                    negh_blk  = surr_blks(BLKNO ,CENTER,RIGHT_EDGE,CENTER,blockID)                    

                    nghproc = -1
                    do p=1,neighProcsCount
                       if (neighProcsList(p) .eq. negh_proc) nghproc = p
                    enddo

                    ! Dof in global numbering:
                    idxb = neighProc_blkList_idx(negh_blk,nghproc)
                    idof_jp1 = ii + (1-1)*NXB + (kk-1)*NXB*NYB + (idxb-1)*dofs_block
                    
!if (gr_meshMe .eq. 0 ) then!.AND. blockID .eq. 1) then
!print*,"Proc",gr_meshMe,"blk",blockID,"ijk",i,j,k,"idxb",idxb,"idof_ip1",idof_jp1
!end if

                    ! Assign to Matrix:
                    nnzglobi         = nnzglobi + 1
                    nnzloci          = nnzloci  + 1
                    nzvali(nnzloci)  = cy                    
                    colindi(nnzloci) = idof_jp1 - 1 !zero-based

!                    !KPD - Incompressible Reference Pressure.
!                    if (lb.eq.1 .AND. i .eq. blkLimits(LOW,IAXIS) .AND. &
!                                      j .eq. blkLimits(LOW,JAXIS) .AND. &
!                                      k .eq. blkLimits(LOW,KAXIS) ) then
!                       nzvali(nnzloci)  = 0.0
!                    end if

                 case(OUTFLOW)

                    ! Add to by:
                    by = by + cy       

                 case(DIRICHLET)

                    ! Add to by:
                    by = by - cy   

                 end select

              else

                    ! Dof in global numbering:
                    nghproc = 1
                    idxb = neighProc_blkList_idx(blockID,nghproc)
                    idof_jp1 = ii + (jj)*NXB + (kk-1)*NXB*NYB + (idxb-1)*dofs_block
                    
!if (gr_meshMe .eq. 0 ) then!.AND. blockID .eq. 1) then
!print*,"Proc",gr_meshMe,"blk",blockID,"ijk",i,j,k,"idxb",idxb,"idof_ip1",idof_jp1
!end if

                    ! Assign to Matrix:
                    nnzglobi         = nnzglobi + 1
                    nnzloci          = nnzloci  + 1
                    nzvali(nnzloci)  = cy                    
                    colindi(nnzloci) = idof_jp1 - 1 !zero-based

!                    !KPD - Incompressible Reference Pressure.
!                    if (lb.eq.1 .AND. i .eq. blkLimits(LOW,IAXIS) .AND. &
!                                      j .eq. blkLimits(LOW,JAXIS) .AND. &
!                                      k .eq. blkLimits(LOW,KAXIS) ) then
!                       nzvali(nnzloci)  = 0.0
!                    end if

              end if


              ! by coefficient
              nghproc = 1
              idxb = neighProc_blkList_idx(blockID,nghproc)
              idof_j = ii + (jj-1)*NXB + (kk-1)*NXB*NYB + (idxb-1)*dofs_block
                    
!print*,"Proc",gr_meshMe,"blk",blockID,"ijk",i,j,k,"idxb",idxb,"idof_j",idof_j

              ! Assign to Matrix:
              nzvali(CONSTANT_ONE)  = nzvali(CONSTANT_ONE) + by                    
              colindi(CONSTANT_ONE) = idof_j - 1 !zero-based              


#if NDIM == MDIM

              ! Discretization coefficients in z:
              !az = 1./(del(KAXIS)**2.0)
              !bz =-2./(del(KAXIS)**2.0)
              !cz = 1./(del(KAXIS)**2.0)
              az = facevarz(MGW8_FACE_VAR,i  ,j  ,k)/(del(KAXIS)**2.0)
              bz =-(facevarz(MGW8_FACE_VAR,i  ,j  ,k)+facevarz(MGW8_FACE_VAR,i  ,j,k+1))/(del(KAXIS)**2.0)
              cz = facevarz(MGW8_FACE_VAR,i  ,j,k+1)/(del(KAXIS)**2.0)



              ! az coefficient:
              if( k .eq. blkLimits(LOW,KAXIS)) then

                 select case(faces(LOW,KAXIS))

                 case(NOT_BOUNDARY,PERIODIC)

                    ! Find the dof of neighbor block in i,j,k-1
                    negh_proc = surr_blks(PROCNO,CENTER,CENTER,LEFT_EDGE,blockID)                    
                    negh_blk  = surr_blks(BLKNO ,CENTER,CENTER,LEFT_EDGE,blockID)                    

                    nghproc = -1
                    do p=1,neighProcsCount
                       if (neighProcsList(p) .eq. negh_proc) nghproc = p
                    enddo

                    ! Dof in global numbering:
                    idxb = neighProc_blkList_idx(negh_blk,nghproc)
                    idof_km1 = ii + (jj-1)*NXB + (NZB-1)*NXB*NYB + (idxb-1)*dofs_block
                    
!if (gr_meshMe .eq. 1 ) then!.AND. blockID .eq. 1) then
!print*,"Proc",gr_meshMe,"blk",blockID,"ijk",i,j,k,"idxb",idxb,"idof_km1",idof_km1
!end if

                    ! Assign to Matrix:
                    nnzglobi         = nnzglobi + 1
                    nnzloci          = nnzloci  + 1
                    nzvali(nnzloci)  = az                    
                    colindi(nnzloci) = idof_km1 - 1 !zero-based

!                    !KPD - Incompressible Reference Pressure.
!                    if (lb.eq.1 .AND. i .eq. blkLimits(LOW,IAXIS) .AND. &
!                                      j .eq. blkLimits(LOW,JAXIS) .AND. &
!                                      k .eq. blkLimits(LOW,KAXIS) ) then
!                       nzvali(nnzloci)  = 0.0
!                    end if
 
                 case(OUTFLOW)

                    ! Add to bz:
                    bz = bz + az                     

                 case(DIRICHLET)

                    ! Add to bz:
                    bz = bz - az   

                 end select


              else

                    ! Dof in global numbering:
                    nghproc = 1
                    idxb = neighProc_blkList_idx(blockID,nghproc)
                    idof_km1 = ii + (jj-1)*NXB + (kk-2)*NXB*NYB + (idxb-1)*dofs_block
                    
!if (gr_meshMe .eq. 1 ) then!.AND. blockID .eq. 1) then
!print*,"Proc",gr_meshMe,"blk",blockID,"ijk",i,j,k,"idxb",idxb,"idof_km1",idof_km1
!end if

                    ! Assign to Matrix:
                    nnzglobi         = nnzglobi + 1
                    nnzloci          = nnzloci  + 1
                    nzvali(nnzloci)  = az                    
                    colindi(nnzloci) = idof_km1 - 1 !zero-based

!                    !KPD - Incompressible Reference Pressure.
!                    if (lb.eq.1 .AND. i .eq. blkLimits(LOW,IAXIS) .AND. &
!                                      j .eq. blkLimits(LOW,JAXIS) .AND. &
!                                      k .eq. blkLimits(LOW,KAXIS) ) then
!                       nzvali(nnzloci)  = 0.0
!                    end if

              endif


              ! cz coefficient:
              if ( k .eq. blkLimits(HIGH,KAXIS)) then

                 select case(faces(HIGH,KAXIS))

                 case(NOT_BOUNDARY,PERIODIC)

                    ! Find the dof of neighbor block in i,j,k+1
                    negh_proc = surr_blks(PROCNO,CENTER,CENTER,RIGHT_EDGE,blockID)                    
                    negh_blk  = surr_blks(BLKNO ,CENTER,CENTER,RIGHT_EDGE,blockID)                    

                    nghproc = -1
                    do p=1,neighProcsCount
                       if (neighProcsList(p) .eq. negh_proc) nghproc = p
                    enddo

                    ! Dof in global numbering:
                    idxb = neighProc_blkList_idx(negh_blk,nghproc)
                    idof_kp1 = ii + (jj-1)*NXB + (1-1)*NXB*NYB + (idxb-1)*dofs_block
                    
!if (gr_meshMe .eq. 0 ) then!.AND. blockID .eq. 1) then
!print*,"Proc",gr_meshMe,"blk",blockID,"ijk",i,j,k,"idxb",idxb,"idof_kp1",idof_kp1
!end if

                    ! Assign to Matrix:
                    nnzglobi         = nnzglobi + 1
                    nnzloci          = nnzloci  + 1
                    nzvali(nnzloci)  = cz                    
                    colindi(nnzloci) = idof_kp1 - 1 !zero-based

!                    !KPD - Incompressible Reference Pressure.
!                    if (lb.eq.1 .AND. i .eq. blkLimits(LOW,IAXIS) .AND. &
!                                      j .eq. blkLimits(LOW,JAXIS) .AND. &
!                                      k .eq. blkLimits(LOW,KAXIS) ) then
!                       nzvali(nnzloci)  = 0.0
!                    end if

                 case(OUTFLOW)

                    ! Add to bz:
                    bz = bz + cz       

                 case(DIRICHLET)

                    ! Add to bz:
                    bz = bz - cz   
              
                 end select

              else

                    ! Dof in global numbering:
                    nghproc = 1
                    idxb = neighProc_blkList_idx(blockID,nghproc)
                    idof_kp1 = ii + (jj-1)*NXB + (kk)*NXB*NYB + (idxb-1)*dofs_block
                    
!if (gr_meshMe .eq. 0 ) then!.AND. blockID .eq. 1) then
!print*,"Proc",gr_meshMe,"blk",blockID,"ijk",i,j,k,"idxb",idxb,"idof_kp1",idof_kp1
!end if

                    ! Assign to Matrix:
                    nnzglobi         = nnzglobi + 1
                    nnzloci          = nnzloci  + 1
                    nzvali(nnzloci)  = cz                    
                    colindi(nnzloci) = idof_kp1 - 1 !zero-based

!                    !KPD - Incompressible Reference Pressure.
!                    if (lb.eq.1 .AND. i .eq. blkLimits(LOW,IAXIS) .AND. &
!                                      j .eq. blkLimits(LOW,JAXIS) .AND. &
!                                      k .eq. blkLimits(LOW,KAXIS) ) then
!                       nzvali(nnzloci)  = 0.0
!                    end if

              end if


              ! bz coefficient
              nghproc = 1
              idxb = neighProc_blkList_idx(blockID,nghproc)
              idof_k = ii + (jj-1)*NXB + (kk-1)*NXB*NYB + (idxb-1)*dofs_block
                    
              ! Assign to Matrix:
              nzvali(CONSTANT_ONE)  = nzvali(CONSTANT_ONE) + bz                    
              colindi(CONSTANT_ONE) = idof_k - 1 !zero-based              


#endif

              
              nzval(nnz_loc+1:nnz_loc+nnzloci)  = nzvali(1:nnzloci)
              colind(nnz_loc+1:nnz_loc+nnzloci) = colindi(1:nnzloci)
              ! Update processor local number of non-zeros:
              nnz_loc = nnz_loc + nnzloci
              

!KPD - For Final Stencil Checking...
!if (gr_meshMe .eq. 0 .AND. blockID .eq. 2 .AND. i.eq.3 .AND. j.eq.3 .AND. k.eq.3) then
!if (gr_meshMe .eq. 0 .AND. blockID .eq. 3 .AND. i.eq.3 .AND. j.eq.3 .AND. k.eq.3) then
!if (gr_meshMe .eq. 0 .AND. blockID .eq. 3 .AND. i.eq.3 .AND. j.eq.10 .AND. k.eq.10) then
!do iii=1,nnzloci
!print*,"Stencil",nrowi,iii,colindi2(iii),nzvali2 (iii)
!end do
!end if

           enddo
        enddo
     enddo

     ! Release Pointers:
     call Grid_releaseBlkPtr(blockID,facevarx ,FACEX)
     call Grid_releaseBlkPtr(blockID,facevary ,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facevarz ,FACEZ)
#endif

  enddo

  rowptr(m_loc+1) = nnz_loc

!  print*,"Final Non-Zero Count on Proc:",gr_meshMe,"is: ",nnz_loc,"first row",fst_row,"# of rows:",m_loc
  write(ind_me,'(I4.4)') gr_meshMe
  open(113,file='./IOData/matrix_'//ind_me//'.res',STATUS='REPLACE')

  do i=1,m_loc
     do j=rowptr(i)+1,rowptr(i+1)
        write(113,*) i+m_loc*gr_meshMe,colind(j)+1,nzval(j)
     enddo
  enddo
  close(113)


  ! Creat the damn beast
  !=====================
  ! m_glob, n_glob = global number of rows
  ! nnz_loc = # of non-zeros on this core
  ! m_loc = # of rows on this core
  ! fst_row = 1st row # on this core (0-based!)
  ! nzval(i) = ith non-zero value on this core
  ! colind(i) = column index of ith nzval on this core (0-based!)
  ! rowptr = non-zero counted entry that starts new row (0-based!)
  !=========================================================================
!if (nnz_loc .gt. 0) then
  call f_dCreate_CompRowLoc_Mat_dist(A, m_glob, n_glob, nnz_loc,           &
                                     m_loc, fst_row, nzval(1:nnz_loc),     &
                                     colind(1:nnz_loc), rowptr(1:m_loc+1), &
                                     SLU_NR_loc, SLU_D, SLU_GE )
!else

!end if
  !=========================================================================

!  if (gr_meshMe .eq. 0) print*,"Leaving gr_sluBuildAdist..." 

  return

end subroutine gr_sluBuildAdist
