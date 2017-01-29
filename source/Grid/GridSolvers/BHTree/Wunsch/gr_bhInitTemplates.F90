!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhInitTemplates
!!
!! NAME
!!
!!  gr_bhInitTemplates
!!
!!
!! SYNOPSIS
!!
!!   gr_bhLeInitTemplates(
!!       integer , intent(in)    :: write_arrs,
!!       integer , intent(inout) :: max_cc_dist,
!!       integer , intent(inout) :: max_cc_ind,
!!       integer , intent(inout) :: max_cn_ind
!!       )
!!
!! DESCRIPTION
!!
!!  Determines interation lists for all templates. There are two 
!!  types of interactions: cell-cell interaction and cell-node interaction.
!!  For each cell in a block, there is a list of cells with which it interacts
!!  (organized into groups with the same distance) and a list of tree-nodes with 
!!  which it interacts.
!!
!!  The subroutine has to be called twice. In the first run, it only determines
!!  maximum numbers of interactions of both types and returns values max_cc_dist,
!!  max_cc_ind and max_cn_ind. Then, arrays gr_bhTreeCn_ind, gr_bhTreeCc_disti, gr_bhTreeCc_count
!!  and gr_bhTreeCc_ind should be allocated. After that, the second call of this 
!!  subroutine fills the arrays with values.
!!
!! ARGUMENTS
!!
!!  write_arrs  - 1: write values into arrays, 0: only determine arrays sizes
!!  max_cc_dist - maximum number of cell-cell distances
!!  max_cc_ind  - maximum number of cell-cell interactions for one distance
!!  max_cn_ind  - maximum number of cell-node interactions
!!
!!***



subroutine gr_bhInitTemplates(write_arrs, max_cc_dist, max_cc_ind, max_cn_ind)

  use gr_bhInterface, ONLY : gr_bhEwald, gr_bhGetTreeSize
  use gr_bhData, ONLY : gr_bhTreeBS, gr_bhTreeLevels, gr_bhTreeLoff, &
    gr_bhTreeCc_count, gr_bhTreeCc_disti, gr_bhTreeCn_ind, gr_bhTreeCc_ind, &
    gr_bhTreeLimAngle2

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer , intent(in)  :: write_arrs
  integer , intent(inout) :: max_cc_dist, max_cc_ind, max_cn_ind
  integer :: i, j, k, l, level, sp
  integer :: ii, jj, kk, cind, cci, cci_last, cni_last
  integer :: temp, xoff, yoff, zoff
  integer :: multi(1:gr_bhTreeLevels), mi(1:gr_bhTreeLevels)
  integer :: fs, fac, pos, dist2, mindist2, ns, nsm1, node_diag2
  integer :: tai, this_dist_i
  integer :: stack(1:gr_bhTreeLevels, 1:(8**(gr_bhTreeLevels+1)-1)/7)
  integer :: cc_dist2(gr_bhTreeBS**3), cc_count(gr_bhTreeBS**3)

  ! ASSUMES BLOCKS ARE CUBES!!! SHOULD BE GENERALIZED

  if (write_arrs .eq. 0) then
    max_cc_dist = 0
    max_cc_ind = 0
    max_cn_ind = 0
  endif
 
  do temp = 1,27
    
    xoff = gr_bhTreeBS * (mod((temp-1), 3) - 1)
    yoff = gr_bhTreeBS * (mod((temp-1), 9) / 3 - 1)
    zoff = gr_bhTreeBS * ((temp-1) / 9 - 1)

    ! FOR EACH CELL OF A TEMPLATE BLOCK WALK THE BLOCK TREE TO FILL INTERACTION LISTS
    do ii = 1, gr_bhTreeBS
      do jj = 1, gr_bhTreeBS
        do kk = 1, gr_bhTreeBS
 
          ! cell index in the interaction list
          cind = (ii-1) + (jj-1)*gr_bhTreeBS + (kk-1)*gr_bhTreeBS**2 + 1
          cci_last = 0 ! counter for cell-cell interaction distances
          cni_last = 0 ! counter for cell-node interaction distances
 
          ! clear the arrays
          if (write_arrs .eq. 1) then
            do i = 1,max_cc_dist
              gr_bhTreeCc_count(cind,i, temp) = -1
            enddo
            do i = 1,max_cn_ind
              gr_bhTreeCn_ind(cind,i, temp) = -1
            enddo
          endif
 
          ! put root node on stack
          do l = 1, gr_bhTreeLevels
            multi(l) = 0
          enddo
          stack(:,1) = multi(:)
          sp = 1
     
          ! walk through the tree
          do
            ! take the node on the bottom of the stack
            multi(:) = stack(:,sp)
            sp = sp - 1
     
            ! determine level of the multi-index
            level = gr_bhTreeLevels
            do l = 1,gr_bhTreeLevels
              if (multi(l) == 0) then
                level = l - 1
                exit
              endif
            enddo
            
            if (level == gr_bhTreeLevels) then
              ! the lowest level => cell-cell interaction
              ! at first convert multi-index to indeces in the block
              i = 1
              j = 1
              k = 1
              do l = 1,gr_bhTreeLevels
                fac = 2**(gr_bhTreeLevels-l)
                i = i + (mod(multi(l)-1,2))      * fac
                j = j + (mod((multi(l)-1)/2,2))  * fac
                k = k + (mod((multi(l)-1)/4,2))  * fac
              enddo
              if ((temp .eq. 14) .and. (i .eq. ii) .and. (j .eq. jj) .and. (k .eq. kk)) goto 10 ! omit the same cell
 
              ! compute position in the gr_bhTreeArray 
              fs = 1
              pos = gr_bhTreeLoff(level)
              do l = level,1,-1
                pos = pos + (multi(l)-1)*fs
                fs = fs * 8
              enddo

              ! compute a distance of the two cells
              dist2 = (i-ii+xoff)*(i-ii+xoff) + (j-jj+yoff)*(j-jj+yoff) + (k-kk+zoff)*(k-kk+zoff)
 
              ! check if this distance already exists in gr_bhTreeCc_disti(cind,:, temp) array
              this_dist_i = -1
              do cci = 1, cci_last
                if (cc_dist2(cci) .eq. dist2) this_dist_i = cci
              enddo
              if (this_dist_i .eq. -1) then
                ! it does not exist => add another distance, set its counter to 1, add the appropriate index 
                cci_last = cci_last + 1
                cc_dist2(cci_last) = dist2
                cc_count(cci_last) = 1
                if (write_arrs .eq. 1) then
                  gr_bhTreeCc_disti(cind, cci_last, temp) = 1./sqrt(real(dist2))
                  gr_bhTreeCc_count(cind, cci_last, temp) = 1
                  gr_bhTreeCc_ind  (cind, cci_last, 1, temp) = pos
                endif
 
              else
                ! it exists => increase gr_bhTreeCc_count and add the index to this cell in gr_bhTreeArray
                cc_count(this_dist_i) = cc_count(this_dist_i) + 1
                if (write_arrs .eq. 1) then
                  gr_bhTreeCc_count(cind, this_dist_i, temp) = gr_bhTreeCc_count(cind, this_dist_i, temp) + 1
                  gr_bhTreeCc_ind  (cind, this_dist_i, gr_bhTreeCc_count(cind, this_dist_i, temp), temp) = pos
                endif
              endif
     
            else
              ! CELL NODE INTERACTION LIST
              ns = 2**(gr_bhTreeLevels-level) ! node size (in cells)
              node_diag2 = 3*ns*ns
              nsm1 = ns - 1
 
              ! find left, back, bottom corner of the node
              i = 1
              j = 1
              k = 1
              do l = 1,level
                fac = 2**(gr_bhTreeLevels-l)
                i = i + (mod(multi(l)-1,2))      * fac
                j = j + (mod((multi(l)-1)/2,2))  * fac
                k = k + (mod((multi(l)-1)/4,2))  * fac
              enddo
              
              ! go through the all 8 corners of the node and find the closest distance
              mindist2 = gr_bhTreeBS**3 + 1
              ! left, back, bottom corner
              dist2 = (i-ii+xoff)*(i-ii+xoff) + (j-jj+yoff)*(j-jj+yoff) + (k-kk+zoff)*(k-kk+zoff)
              if (dist2 .lt. mindist2) mindist2 = dist2
              ! right, back, bottom corner
              dist2 = ((i+nsm1)-ii+xoff)*((i+nsm1)-ii+xoff) + (j-jj+yoff)*(j-jj+yoff) + (k-kk+zoff)*(k-kk+zoff)
              if (dist2 .lt. mindist2) mindist2 = dist2
              ! left, front, bottom corner
              dist2 = (i-ii+xoff)*(i-ii+xoff) + ((j+nsm1)-jj+yoff)*((j+nsm1)-jj+yoff) + (k-kk+zoff)*(k-kk+zoff)
              if (dist2 .lt. mindist2) mindist2 = dist2
              ! right, front, bottom corner
              dist2 = ((i+nsm1)-ii+xoff)*((i+nsm1)-ii+xoff) + ((j+nsm1)-jj+yoff)*((j+nsm1)-jj+yoff) + (k-kk+zoff)*(k-kk+zoff)
              if (dist2 .lt. mindist2) mindist2 = dist2
              
              ! left, back, top corner
              dist2 = (i-ii+xoff)*(i-ii+xoff) + (j-jj+yoff)*(j-jj+yoff) + ((k+nsm1)-kk+zoff)*((k+nsm1)-kk+zoff)
              if (dist2 .lt. mindist2) mindist2 = dist2
              ! right, back, top corner
              dist2 = ((i+nsm1)-ii+xoff)*((i+nsm1)-ii+xoff) + (j-jj+yoff)*(j-jj+yoff) + ((k+nsm1)-kk+zoff)*((k+nsm1)-kk+zoff)
              if (dist2 .lt. mindist2) mindist2 = dist2
              ! left, front, top corner
              dist2 = (i-ii+xoff)*(i-ii+xoff) + ((j+nsm1)-jj+yoff)*((j+nsm1)-jj+yoff) + ((k+nsm1)-kk+zoff)*((k+nsm1)-kk+zoff)
              if (dist2 .lt. mindist2) mindist2 = dist2
              ! right, front, top corner
              dist2 = ((i+nsm1)-ii+xoff)*((i+nsm1)-ii+xoff) + ((j+nsm1)-jj+yoff)*((j+nsm1)-jj+yoff) &
                                                            + ((k+nsm1)-kk+zoff)*((k+nsm1)-kk+zoff)
              if (dist2 .lt. mindist2) mindist2 = dist2
 

              ! test size/dist < angle
              if (node_diag2 .lt. gr_bhTreeLimAngle2*mindist2) then
                cni_last = cni_last + 1
 
                if (write_arrs .eq. 1) then
                  ! find the index in the tree array
                  fs = 4
                  tai = gr_bhTreeLoff(level)
                  do l = level,1,-1
                    tai = tai + (multi(l)-1)*fs
                    fs = fs * 8
                  enddo
                  gr_bhTreeCn_ind(cind, cni_last, temp) = tai
                endif
 
              else
                ! put all children on the stack
                mi(:) = multi(:)
                do i = 1,8
                  mi(level+1) = i
                  sp = sp + 1
                  stack(:,sp) = mi(:)
                enddo
              endif
            endif
 
10          continue ! end 
            ! stack is empty - exiting
            if (sp == 0) exit
          enddo

          if (cci_last .gt. max_cc_dist) max_cc_dist = cci_last
          if (cni_last .gt. max_cn_ind ) max_cn_ind  = cni_last
          do i = 1, cci_last
            if (cc_count(i) .gt. max_cc_ind) then
              max_cc_ind = cc_count(i)
            endif
          enddo
 
        enddo
      enddo
    enddo
    if (write_arrs .eq. 0) then
      max_cc_dist = max_cc_dist + 1
      max_cn_ind  = max_cn_ind  + 1
    endif
  enddo

  return
end


