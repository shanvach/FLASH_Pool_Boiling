!!****if* source/Simulation/SimulationMain/magnetoHD/RTmhd/Grid_bcApplyToRegionSpecialized
!!
!!
!! NAME
!!  Grid_bcApplyToRegionSpecialized
!!
!! SYNOPSIS
!!
!!  Grid_bcApplyToRegionSpecialized(integer(IN)  :: bcType,
!!                                 integer(IN)  :: gridDataStruct,
!!                                 integer(IN)  :: guard,
!!                                 integer(IN)  :: axis,
!!                                 integer(IN)  :: face,
!!                                 real(INOUT)  :: regionData(:,:,:,:),
!!                                 integer(IN)  :: regionSize(:),
!!                                 logical(IN)  :: mask(:),
!!                                 logical(OUT) :: applied,
!!                                 integer(IN)  :: blockHandle,secondDir,thirdDir,
!!                                 integer(IN)  :: endPoints(LOW:HIGH,MDIM),
!!                                 integer(IN)  :: blkLimitsGC(LOW:HIGH,MDIM),
!!                        OPTIONAL,integer(IN)  :: idest )
!!
!!                    
!!  
!! DESCRIPTION 
!!  
!!  Applies the boundary conditions to the specified data structure.
!!  The routine is handed a region that has been extracted from the
!!  data structure, on which it should apply the boundary conditions. 
!!  The direction along which the BC are to be applied is always the first
!!  dimension in the given region, and the last dimension contains the
!!  the variables in the data structure. The middle two dimension contain
!!  the size of the region along the two dimensions of the physical grid
!!  that are not having the BC applied.
!!
!!  This routine applies the boundary conditions on a given face (lowerface
!!  or upperface) along a given axis, by using and setting values
!!  for all variables in the gridDataStruct that are not masked out. The 
!!  argument "mask" has the information about the masked variables.
!! 
!!     If (face=LOW)  
!!       regionData(1:guard,:,:,masked(variables) =  boundary values
!!     If (face=HIGH) 
!!       regionData(regionSize(BC_DIR)-guard+1:regionSize(BC_DIR),:,:,masked(variables) =  boundary values
!!
!!
!!   This interface serves two purposes. One it provides a mechanism by which users can
!!   supply their own custom boundary conditions, and two it allows FLASH to implement
!!   more complex boundary conditions such as those with hydrostatic equilibrium in
!!   isolation from the machinery of setting up the region etc. This interface has
!!   two extra arguments over Grid_bcApplyToRegion. One is the blockHandle, which allows
!!   access to the coordinates information, and second one is logical variable "applied".
!!   This routine is always called first when handling boundary conditions, and if it
!!   finds a match with one of bcTypes that it implements, it sets "applied" to .true., otherwise
!!   "applies" is set to .false. If applied is false, then Grid_bcApplyToRegion is called, and
!!   that is the only routine out of these two that traps the error if none of the boundary
!!   types are matched in either routine.
!!
!! ARGUMENTS 
!!
!!  bcType - the type of boundary conditions being applied
!!  gridDataStruct - the Grid dataStructure
!!  guard -    number of guardcells 
!!  axis  - the dimension along which to apply boundary conditions,
!!          can take values of IAXIS, JAXIS and KAXIS
!!  face    -  can take values LOW and HIGH, defined in constants.h
!!             to indicate whether to apply boundary on lowerface or 
!!             upperface
!!  regionData     : the extracted region from a block of permanent storage of the 
!!                   specified data structure. Its size is given by regionSize
!!  regionSize     : regionSize(BC_DIR) contains the size of the each row and
!!                   regionSize(SECOND_DIR) contains the number of along the second
!!                   direction, and regionSize(THIRD_DIR) has the number of rows
!!                   along the third direction. regionSize(GRID_DATASTRUCT) contains the
!!                   number of variables in the data structure
!!  mask - if present applies boundary conditions to only selected variables
!!  blockHandle - the local block number
!!  applied - is set true if this routine matched a bcType, otherwise it is set to false
!!
!!  idest - Only meaningful with PARAMESH 3 or later.  The argument indicates which slot
!!          in its one-block storage space buffers ("data_1blk.fh") PARAMESH is in the
!!          process of filling.
!!          The following applies when guard cells are filled as part of regular
!!          Grid_fillGuardCells processing (or, in NO_PERMANENT_GUARDCELLS mode,
!!          in order to satisfy a Grid_getBlkPtr request): The value is 1 if guard cells
!!          are being filled in the buffer slot in UNK1 and/or FACEVAR{X,Y,Z}1 or WORK1
!!          that will end up being copied to permanent block data storage (UNK and/or
!!          FACEVAR{X,Y,Z} or WORK, respectively) and/or returned to the user.
!!          The value is 2 if guard cells are being filled in the alternate slot in
!!          the course of assembling data to serve as input for coarse-to-fine
!!          interpolation.
!!          When guard cells are being filled in order to provide input data for
!!          coarse-to-fine interpolation as part of amr_prolong processing (which
!!          is what happens when Grid_updateRefinement is called for an AMR Grid),
!!          the value is always 1.
!!
!!          In other words, you nearly always want to ignore this optional
!!          argument.  This implementation does just that.
!!
!!  NOTES 
!!
!!            This routine is common to all the mesh packages supported
!!            The mesh packages extract the small vectors relevant to
!!            boundary conditions calculations from their Grid data 
!!            structures. 
!!            
!!
!!***


subroutine Grid_bcApplyToRegionSpecialized(bcType,gridDataStruct,&
     guard,axis,face,regionData,regionSize,mask,&
     applied,blockHandle,secondDir,thirdDir,endPoints,blkLimitsGC, idest)

#include "constants.h"
#include "Flash.h"

  use Grid_data,        ONLY : gr_meshMe,gr_domainBC
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface,   ONLY : Grid_getBlkPtr,Grid_releaseBlkPtr
  use Gravity_data,     ONLY : grv_const
  !use gr_bcInterface,   ONLY : gr_bcMapBcType

  implicit none

  integer, intent(IN) :: bcType,axis,face,guard,gridDataStruct
  integer,intent(IN) :: secondDir,thirdDir
  integer,dimension(REGION_DIM),intent(IN) :: regionSize
  real,dimension(regionSize(BC_DIR),&
       regionSize(SECOND_DIR),&
       regionSize(THIRD_DIR),&
       regionSize(STRUCTSIZE)),intent(INOUT)::regionData
  logical,intent(IN),dimension(regionSize(STRUCTSIZE)):: mask
  integer,intent(IN) :: blockHandle
  integer,intent(IN),dimension(LOW:HIGH,MDIM) :: endPoints, blkLimitsGC
  logical, intent(OUT) :: applied
  integer,intent(IN),OPTIONAL:: idest

  logical,dimension(regionSize(STRUCTSIZE)) :: localMask
  integer :: i,j, k,ivar,ind,je,ke,n,varCount
  logical :: isFace
  real    :: sign

  integer :: sizeX,sizeY,istat,xOffset,zOffset,ig
  real, allocatable, dimension(:) :: yCoord


  !! Set this true if user boundary condition is needed
  applied=.true.


  je=regionSize(SECOND_DIR)
  ke=regionSize(THIRD_DIR)
  varCount=regionSize(STRUCTSIZE)
  isFace = (gridDataStruct==FACEX).and.(axis==IAXIS)
  isFace = isFace.or.((gridDataStruct==FACEY).and.(axis==JAXIS))
  isFace = isFace.or.((gridDataStruct==FACEZ).and.(axis==KAXIS))


  call gr_bcGetCoords_internal

  xOffset = endPoints(LOW,secondDir) - 1
  zOffset = endPoints(LOW,thirdDir) - 1


  if (bcType .ne. USER_DEFINED) then
     do ivar = 1,varCount

        if(mask(ivar)) then
           !call gr_bcMapBcType(bcTypeActual,bcType,ivar,gridDataStruct,axis,face,idest)
           sign = 1.
           if (gridDataStruct==CENTER) then
#ifdef VELX_VAR
              if ((axis==IAXIS).and.(ivar==VELX_VAR))sign=-1.
#endif
#ifdef VELY_VAR
              if((axis==JAXIS).and.(ivar==VELY_VAR))sign=-1.
#endif
#ifdef VELZ_VAR
              if((axis==KAXIS).and.(ivar==VELZ_VAR))sign=-1.
#endif
           end if


           if(face==LOW) then ! apply BC on the bottom boundary at y=0

              select case (bcType)
              case(REFLECTING)
                 k = 2*guard+1
                 if(isFace)k=k+1
                 do i = 1,guard
                    regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)*sign
                 end do


              case(DIODE)
                 do i = 1,guard
                    regionData(i,1:je,1:ke,ivar)= regionData(guard+1,1:je,1:ke,ivar)
                    if (sign == -1) then
                       do n=1,ke
                          do j=1,je
                             regionData(i,j,n,ivar) = min(regionData(i,j,n,ivar),0.0)
                          end do
                       end do
                    end if
                 end do


              case(DIRICHLET)
                 do i = 1,guard
                    regionData(guard+1-i,1:je,1:ke,ivar)= (1-2*i)*regionData(guard+1,1:je,1:ke,ivar)
                 end do


              case(GRIDBC_MG_EXTRAPOLATE)
                 do i = 1,guard
                    regionData(guard+1-i,1:je,1:ke,ivar)= (1+i)*regionData(guard+1,1:je,1:ke,ivar) &
                         - i*regionData(guard+2,1:je,1:ke,ivar)
                 end do


              case(OUTFLOW,HYDROSTATIC_F2_NVOUT,HYDROSTATIC_F2_NVDIODE,HYDROSTATIC_F2_NVREFL)
                 do i = 1,guard
                    regionData(i,1:je,1:ke,ivar)= regionData(guard+1,1:je,1:ke,ivar)
                 end do


              case default
                 print*,'boundary is',bcType
                 call Driver_abortFlash("unsupported boundary condition on Lower Face")
              end select


           else !apply BC on the top boundary at y=ymax

              select case (bcType)
              case(REFLECTING)
                 k = 2*guard+1
                 if(isFace)k=k+1
                 do i = 1,guard
                    regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)*sign
                 end do


              case(DIODE)
                 k=guard
                 if(isFace)k=k+1
                 do i = 1,guard
                    regionData(k+i,1:je,1:ke,ivar)= regionData(k,1:je,1:ke,ivar)
                    if (sign == -1) then
                       do n = 1,ke
                          do j = 1,je
                             regionData(k+i,j,n,ivar) = max(regionData(k+i,j,n,ivar),0.0)
                          end do
                       end do
                    end if
                 end do


              case(DIRICHLET)
                 k=guard
                 do i = 1,guard
                    regionData(k+i,1:je,1:ke,ivar)= (1-2*i)*regionData(k,1:je,1:ke,ivar)
                 end do


              case(GRIDBC_MG_EXTRAPOLATE)
                 k=guard
                 do i = 1,guard
                    regionData(k+i,1:je,1:ke,ivar)= (1+i)*regionData(k,1:je,1:ke,ivar) &
                         - i*regionData(k-1,1:je,1:ke,ivar)
                 end do


              case(OUTFLOW,HYDROSTATIC_F2_NVOUT,HYDROSTATIC_F2_NVDIODE,HYDROSTATIC_F2_NVREFL)
                 k=guard
                 if(isFace)k=k+1
                 do i = 1,guard
                    regionData(k+i,1:je,1:ke,ivar)= regionData(k,1:je,1:ke,ivar)
                 end do


              case default
                 print*,'boundary is',bcType
                 call Driver_abortFlash("unsupported boundary condition on Upper Face")
              end select
           end if
        end if
     end do !end of do i=1,varCount


  elseif (bcType == USER_DEFINED) then

     if(face==LOW) then
        ig=guard
        do i=1,guard
           !! First copy everything as outflow BC
           k=guard
           if(isFace)k=k+1
           regionData(k+i,1:je,1:ke,1:varCount)= regionData(k,1:je,1:ke,1:varCount)
        enddo

        do i = 1,guard !(this would be indices for interior cells)
           k=2*guard+1
           !! Modify pressure and energy using the hydrostatic BC
           if (gridDataStruct==CENTER) then
              !! PRESSURE
              regionData(k-i,1:je,1:ke,PRES_VAR) = regionData(ig,1:je,1:ke,PRES_VAR)+&
                   grv_const*regionData(k-i,1:je,1:ke,DENS_VAR)*abs(yCoord(k-i)-yCoord(ig))

              !! ENERGY
              regionData(k-i,1:je,1:ke,ENER_VAR) = regionData(k-i,1:je,1:ke,PRES_VAR)&
                   /(regionData(k-i,1:je,1:ke,DENS_VAR)*(regionData(k-i,1:je,1:ke,GAME_VAR)-1.))&
                   +0.5*( regionData(k-i,1:je,1:ke,VELX_VAR)**2 &
                         +regionData(k-i,1:je,1:ke,VELY_VAR)**2 &
                         +regionData(k-i,1:je,1:ke,VELZ_VAR)**2 )
           endif
        end do

     else ! face==HIGH
        ig=guard
        do i=1,guard
           !! First copy everything as outflow BC
           k=guard
           if(isFace)k=k+1
           regionData(k+i,1:je,1:ke,1:varCount)= regionData(k,1:je,1:ke,1:varCount)
        enddo

        do i = 1,guard !(this would be indices for interior cells)
           k=2*guard+1
           !! Modify pressure and energy using the hydrostatic BC
           if (gridDataStruct==CENTER) then
              !! PRESSURE
              regionData(k-i,1:je,1:ke,PRES_VAR) = regionData(ig,1:je,1:ke,PRES_VAR)+&
                   grv_const*regionData(k-i,1:je,1:ke,DENS_VAR)*abs(yCoord(k-i)-yCoord(ig))

              !! ENERGY
              regionData(k-i,1:je,1:ke,ENER_VAR) = regionData(k-i,1:je,1:ke,PRES_VAR)&
                   /(regionData(k-i,1:je,1:ke,DENS_VAR)*(regionData(k-i,1:je,1:ke,GAME_VAR)-1.))&
                   +0.5*( regionData(k-i,1:je,1:ke,VELX_VAR)**2 &
                         +regionData(k-i,1:je,1:ke,VELY_VAR)**2 &
                         +regionData(k-i,1:je,1:ke,VELZ_VAR)**2 )
           endif
        end do
     end if !end if (face==LOW)
  end if !end if (bcType == USER_DEFINED) then





  !! deallocate
  deallocate(yCoord)

  return

contains
  subroutine gr_bcGetCoords_internal
    use Grid_data, ONLY : gr_meshMe
    integer :: sizeGC

    !! This is for the first coordinate that is along the direction BC is applied to
    sizeGC = blkLimitsGC(HIGH,axis)

    !! allocate
    allocate(yCoord(sizeGC))
    if (isFace) then
       call gr_extendedGetCellCoords(axis, blockHandle, gr_meshMe, FACES, .true., yCoord, sizeGC)
    else
       call gr_extendedGetCellCoords(axis, blockHandle, gr_meshMe, CENTER, .true., yCoord, sizeGC)
    end if

!!$    !! This is for the second coordinate that is transverse the direction BC is applied to
!!$    sizeGC = blkLimitsGC(HIGH,secondDir)
!!$    allocate(xCoord(sizeGC))
!!$#if NDIM > 1
!!$    if ((gridDataStruct==FACEX .AND. secondDir==IAXIS) .OR. &
!!$         (gridDataStruct==FACEY .AND. secondDir==JAXIS) .OR. &
!!$         (gridDataStruct==FACEZ .AND. secondDir==KAXIS)) then
!!$       call gr_extendedGetCellCoords(secondDir, blockHandle, gr_meshMe, FACES, .true., xCoord, sizeGC)
!!$    else
!!$       call gr_extendedGetCellCoords(secondDir, blockHandle, gr_meshMe, CENTER, .true., xCoord, sizeGC)
!!$    end if
!!$#endif
!!$    sizeGC = blkLimitsGC(HIGH,thirdDir)
!!$    allocate(thirdCoord(sizeGC))
!!$#if NDIM > 2
!!$    if ((gridDataStruct==FACEX .AND. thirdDir==IAXIS) .OR. &
!!$         (gridDataStruct==FACEY .AND. thirdDir==JAXIS) .OR. &
!!$         (gridDataStruct==FACEZ .AND. thirdDir==KAXIS)) then
!!$       call gr_extendedGetCellCoords(thirdDir, blockHandle, gr_meshMe, FACES, .true., thirdCoord, sizeGC)
!!$    else
!!$       call gr_extendedGetCellCoords(thirdDir, blockHandle, gr_meshMe, CENTER, .true., thirdCoord, sizeGC)
!!$    end if
!!$#endif

  end subroutine gr_bcGetCoords_internal



end subroutine Grid_bcApplyToRegionSpecialized
