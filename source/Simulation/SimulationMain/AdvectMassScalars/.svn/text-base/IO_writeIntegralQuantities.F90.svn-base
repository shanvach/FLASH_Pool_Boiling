!!****if* source/Simulation/SimulationMain/AdvectMassScalars/IO_writeIntegralQuantities
!!
!!
!!  NAME
!!    IO_writeIntegralQuantities
!!
!!  SYNOPSIS
!!    call IO_writeIntegralQuantities() 
!!                                    integer(in) :: isFirst,
!!                                    real(in)    :: simTime)
!!
!!  DESCRIPTION
!!
!!   Compute the values of integral quantities (eg. total energy)
!!   and write them to an ASCII file.  If this is the initial step,
!!   create the file and write a header to it before writing the data.
!!
!!   Presently, this supports 1, 2, and 3-d Cartesian geometry and 2-d
!!   cylindrical geometry (r,z).  More geometries can be added by
!!   modifying the volume of each zone (dvol).
!!
!!   Users should modify this routine if they want to store any
!!   quantities other than default values in the flash.dat.  Make sure
!!   to modify the nGlobalSum parameter to match the number of
!!   quantities written.  Also make sure to modify the header to match
!!   the names of quantities with those calculated in the lsum and
!!   gsum arrays.
!!  
!!  ARGUMENTS
!!    
!!   io_globalMe  - current processor number
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!   simTime - simulation time
!!
!!
!!***

subroutine IO_writeIntegralQuantities ( isFirst, simTime)

  use IO_data, ONLY : io_restart, io_statsFileName
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_getDeltas, &
    Grid_releaseBlkPtr

   use IO_data, ONLY : io_globalMe
  implicit none

#include "mpif.h"
#include "constants.h"
#include "Flash.h"
  
  
  real, intent(in) :: simTime

  integer, intent(in) :: isFirst

  integer :: lb, count
  
  integer :: funit = 99
  integer :: error
  
  character (len=MAX_STRING_LENGTH), save :: fname 
  
  integer :: blockList(MAXBLOCKS)

  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)

  integer, parameter ::  nGlobalSum = 12          ! Number of globally-summed quantities
  real :: gsum(nGlobalSum) !Global summed quantities
  real :: lsum(nGlobalSum) !Global summed quantities

  integer :: i, j, k
  real :: dvol, del(MDIM)
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  integer :: point(MDIM)

  ! Sum quantities over all locally held leaf-node blocks.
  gsum  = 0.
  lsum = 0.
  
  call Grid_getListOfBlocks(LEAF, blockList, count)
  
  do lb = 1, count
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)



     !DEV: only in place for cartesian coords
     !getting the cell volume
     call Grid_getDeltas(blocklist(lb), del)

     if(NDIM == 1) then
        dvol = del(IAXIS)
     else if(NDIM == 2) then
        dvol = del(IAXIS) * del(JAXIS)
     else
        dvol = del(IAXIS) * del(JAXIS) * del(KAXIS)
     end if



     ! Sum contributions from the indicated range of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k
              
!! DEV: For non-cartesian coordinates, we could used this here instead of the
!! explicit computation of dvol above:             
!               call Grid_getSingleCellVol(blockList(i), EXTERIOR, point, dvol)
     
              ! mass   
#ifdef DENS_VAR
              lsum(1) = lsum(1) + solnData(DENS_VAR,i,j,k)*dvol 
#endif           


#ifdef DENS_VAR
#ifdef VELX_VAR      
              ! momentum
              lsum(2) = lsum(2) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELX_VAR,i,j,k)*dvol 
           
#endif
#ifdef VELY_VAR      

              lsum(3) = lsum(3) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELY_VAR,i,j,k)*dvol
           
#endif
#ifdef VELZ_VAR      
              lsum(4) = lsum(4) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELZ_VAR,i,j,k)*dvol
#endif

              ! total energy
#ifdef ENER_VAR
              lsum(5) = lsum(5) + solnData(ENER_VAR,i,j,k) * & 
                   &                                solnData(DENS_VAR,i,j,k)*dvol 
#endif

           
#ifdef VELX_VAR      
#ifdef VELY_VAR      
#ifdef VELZ_VAR      
              ! kinetic energy
              lsum(6) = lsum(6) + 0.5*solnData(DENS_VAR,i,j,k) * & 
                   &                             (solnData(VELX_VAR,i,j,k)**2+ & 
                   &                              solnData(VELY_VAR,i,j,k)**2+ & 
                   &                              solnData(VELZ_VAR,i,j,k)**2)*dvol           

#endif
#endif
#endif


#ifdef EINT_VAR
              ! internal energy
              lsum(7) = lsum(7) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(EINT_VAR,i,j,k)*dvol
#endif

#ifdef MS1_MSCALAR
              lsum(8) = lsum(8) + solnData(DENS_VAR,i,j,k) * dvol * solnData(MS1_MSCALAR,i,j,k)
#endif

#ifdef MS2_MSCALAR
              lsum(9) = lsum(9) + solnData(DENS_VAR,i,j,k) * dvol * solnData(MS2_MSCALAR,i,j,k)
#endif

#ifdef MS3_MSCALAR
              lsum(10) = lsum(10) + solnData(DENS_VAR,i,j,k) * dvol * solnData(MS3_MSCALAR,i,j,k)
#endif
           
#ifdef MS4_MSCALAR
              lsum(11) = lsum(11) + solnData(DENS_VAR,i,j,k) * dvol * solnData(MS4_MSCALAR,i,j,k)
#endif
           
#ifdef MS5_MSCALAR
              lsum(12) = lsum(12) + solnData(DENS_VAR,i,j,k) * dvol * solnData(MS5_MSCALAR,i,j,k)
#endif
#endif
           
           enddo
        enddo
     enddo
     call Grid_releaseBlkPtr(blockList(lb), solnData)

  enddo
  

  
  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Reduce (lsum, gsum, nGlobalSum, MPI_Double_Precision, MPI_Sum, & 
       &                MASTER_PE, MPI_Comm_World, error)
  

  if (io_globalMe  == MASTER_PE) then
     
     ! create the file from scratch if it is a not a restart simulation, 
     ! otherwise append to the end of the file
     if (isfirst == 0) then
        open (funit, file=trim(io_statsFileName), position='APPEND')
     else 
        if (.NOT. io_restart) then
           open (funit, file=trim(io_statsFileName)) 
           write (funit, 10)               &
                '#time                     ', &
                'mass                      ', &
                'x-momentum                ', &
                'y-momentum                ', & 
                'z-momentum                ', &
                'E_total                   ', &
                'E_kinetic                 ', &
                'E_internal                ', & 
                'mass (ms 1)               ', &
                'mass (ms 2)               ', &
                'mass (ms 3)               ', &
                'mass (ms 4)               ', &
                'mass (ms 5)               '
           
           
10         format (2x,50(a25, :, 1X))
           
        else
           open (funit, file=trim(io_statsFileName), position='APPEND')
           write (funit, 11) 
11         format('# simulation restarted')
        endif
     endif
     
     write (funit, 12) simtime, gsum      ! Write the global sums to the file.
12   format (1x, 50(es25.18, :, 1x))
 
     close (funit)          ! Close the file.
     
  endif
  
  call MPI_Barrier (MPI_Comm_World, error)
  
  !=============================================================================
  
  return
end subroutine IO_writeIntegralQuantities



