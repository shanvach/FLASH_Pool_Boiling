! Output the information of the input arraies.
! For debug
! Shizhao Wang
! Jan 12, 2015
! =============================================================================
#include "Flash.h"
#include "constants.h"

      subroutine check_array(nxa,nya,nza,nxb,nyb,nzb,nxc,nyc,nzc,a,b,c,id,note)
      implicit none
      integer :: nxa, nya, nza
      integer :: nxb, nyb, nzb
      integer :: nxc, nyc, nzc
      real :: a(nxa,nya,nza), b(nxb,nyb,nzb), c(nxc,nyc,nzc)
      integer :: id
      character(len=*) :: note

      write(10000+id,*) note,', ', id
      write(10000+id,*) 'nxa, nya, nza:', nxa, nya, nza
      write(10000+id,*) 'suma, sumb, id:', sum(a(GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)), & 
                                           sum(b(GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI)), &
                                           sum(c(GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI+1)), id
      write(10000+id,*) 'maxa, maxb, id:', maxval(a(GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)), & 
                                           maxval(b(GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI)), &
                                           maxval(c(GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI+1)), id
      write(10000+id,*) 'mina, minb, id:', minval(a(GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)), & 
                                           minval(b(GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI)), &
                                           minval(c(GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI)), id
      
      return
      end subroutine check_array
! =============================================================================
! =============================================================================


