!!****if* source/Simulation/SimulationMain/Nuc2Grid/Grid_bcApplyToRegion
!!
!! NAME
!!  Grid_bcApplyToRegion
!!
!! SYNOPSIS
!!
!!  call Grid_bcApplyToRegion(integer(IN)  :: bcType,
!!                            integer(IN)  :: gridDataStruct,
!!                            integer(IN)  :: guard,
!!                            integer(IN)  :: axis,
!!                            integer(IN)  :: face,
!!                            real(INOUT)  :: regionData(:,:,:,:),
!!                            integer(IN)  :: regionSize(:),
!!                            logical(IN)  :: mask(:),
!!                            logical(OUT) :: applied,
!!                            integer(IN)  :: blockHandle,
!!                            integer(IN)  :: secondDir,
!!                            integer(IN)  :: thirdDir,
!!                            integer(IN)  :: endPoints(LOW:HIGH,MDIM),
!!                            integer(IN)  :: blkLimitsGC(LOW:HIGH,MDIM),
!!                   OPTIONAL,integer(IN)  :: idest)
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
!!  This version adds handling for USER_DEFINED boundary conditions which
!!  is somewhat like REFLECTING without the special handling for normal velocity
!!  components, with additional handling for species mass fractions.
!!
!!  This routine applies the boundary conditions on a given face (lowerface
!!  or upperface) along a given axis, by using and setting values
!!  for all variables in the gridDataStruct that are not masked out. The 
!!  argument "mask" has the information about the masked variables.
!! 
!!     If (face=LOW)  
!!       regionData(1:guard,:,:,masked(variables) =  boundary values
!!     If (face=HIGH) 
!!       regionData(regionSize(BC_DIR)-guard+1:regionSize(BCDIR),:,:,masked(variables) =  boundary values
!!
!!  One reason why information about direction and variable is
!!  included in this interface is because Velocities need to be
!!  treated specially for REFLECTING boundary conditions. if
!!  axis=IAXIS, then the variable VELX_VAR is treated differently,
!!  same with VELY_VAR if axis=JAXIS and VELZ_VAR if
!!  axis=KAXIS. All supported mesh packages extract the vector passed
!!  in through the argument "dataRow" from the appropriated blocks,
!!  and send it to this routine for boundary calculation. The
!!  PERIODIC boundary is calculated by default when the blocks are
!!  exchanging data with each other.
!!  This routine currently passes handling of all other boundary condition
!!  types on to Grid_applyBCEdge, which is called for each of the variables
!!  in turn.
!!
!!  This routine support only simple boundary conditions that are applied strictly
!!  directionally, and have no need for other grid information such as the coordinates etc.
!!  Currently supported boundary conditions are "OUTFLOW", "REFLECTING" and "DIODE".
!!  The "PERIODIC" boundary conditions are automatically applied in the process of filling
!!  the guard cells all over the domain, and therefore do not need to call this routine
!!
!!  If the user wishes to apply different boundary conditions, they can either
!!  use the interface Grid_bcApplyToRegionSpecialized, or make a copy of this
!!  routine in their simulation directory and customize it.
!!
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
!!  applied - is set true if this routine has handled the given bcType, otherwise it is 
!!            set to false.
!!
!!  blockHandle - Handle for the block for which guardcells are to be filled.
!!              In grid implementations other than Paramesh 4, this is always
!!              a local blockID.
!!
!!              With Paramesh 4:
!!              This may be a block actually residing on the local processor,
!!              or the handle may refer to a block that belong to a remote processor
!!              but for which cached information is currently available locally.
!!              The two cases can be distinguished by checking whether 
!!              (blockHandle .LE. lnblocks): this is true only for blocks that
!!              reside on the executing processor.
!!              The block ID is available for passing on to some handlers for 
!!              boundary conditions that may need it, ignored in the default 
!!              implementation.
!!
!!  secondDir,thirdDir -   Second and third coordinate directions.
!!                         These are the transverse directions perpendicular to
!!                         the sweep direction.
!!                         This is not needed for simple boundary condition types
!!                         such as REFLECTIVE or OUTFLOW, It is provided for
!!                         convenience so that more complex boundary condition
!!                         can make use of it.
!!                         The values are currently fully determined by the sweep
!!                         direction bcDir as follows:
!!                          bcDir   |    secondDir       thirdDir
!!                          ------------------------------------------
!!                          IAXIS   |    JAXIS             KAXIS
!!                          JAXIS   |    IAXIS             KAXIS
!!                          KAXIS   |    IAXIS             JAXIS
!!
!!  endPoints - starting and endpoints of the region of interest.
!!              See also NOTE (1) below.
!!
!!  blkLimitsGC - the starting and endpoint of the whole block including
!!                the guard cells, as returned by Grid_getBlkIndexLimits.
!!              See also NOTE (1) below.
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
!!          argument.  As of FLASH 3.1, it is only used internally within the
!!          Grid unit by a Multigrid GridSolver implementation.
!!
!!  NOTES 
!!
!!            This routine is common to all the mesh packages supported
!!            The mesh packages extract the small vectors relevant to
!!            boundary conditions calculations from their Grid data 
!!            structures. If users wish to apply a different boundary condition, 
!!            they should look at routine Grid_bcApplyToRegionSpecialzed.
!!            They can replace this routine, the implentations included apply
!!            some hydrostatic boundary conditions.
!!            
!!
!!***

subroutine Grid_bcApplyToRegion(bcType,gridDataStruct,&
          guard,axis,face,regionData,regionSize,mask,applied,&
     blockHandle,secondDir,thirdDir,endPoints,blkLimitsGC, idest)

#include "constants.h"
#include "Flash.h"

  use Driver_interface, ONLY : Driver_abortFlash
  use gr_bcInterface, ONLY : gr_bcMapBcType

  implicit none
  
  integer, intent(IN) :: bcType,axis,face,guard,gridDataStruct
  integer,dimension(REGION_DIM),intent(IN) :: regionSize
  real,dimension(regionSize(BC_DIR),&
       regionSize(SECOND_DIR),&
       regionSize(THIRD_DIR),&
       regionSize(STRUCTSIZE)),intent(INOUT)::regionData
  logical,intent(IN),dimension(regionSize(STRUCTSIZE)):: mask
  logical, intent(OUT) :: applied
  integer,intent(IN) :: blockHandle
  integer,intent(IN) :: secondDir,thirdDir
  integer,intent(IN),dimension(LOW:HIGH,MDIM) :: endPoints, blkLimitsGC
  integer,intent(IN),OPTIONAL:: idest

  integer :: i,j, k,ivar,je,ke,n,varCount,bcTypeActual
  logical :: isFace
  integer    :: sign
  real    :: fixVal

  select case (bcType)
  case(REFLECTING,USER_DEFINED,OUTFLOW,DIODE)
     applied = .TRUE.           !will handle these types of BCs below
  case(HYDROSTATIC_F2_NVOUT,HYDROSTATIC_F2_NVDIODE,HYDROSTATIC_F2_NVREFL, &
       HYDROSTATIC_NVOUT,HYDROSTATIC_NVDIODE,HYDROSTATIC_NVREFL)
     if (gridDataStruct==CENTER) then
        applied = .FALSE.       !should be picked up by Flash2HSE implementation
                                !or Flash3HSE implementation if included.
        return                  !RETURN immediately!
     else
        applied = .TRUE.           !will handle these types below (like OUTFLOW)
     end if
  case default
     applied = .FALSE.
     return                     !RETURN immediately!
  end select


  je=regionSize(SECOND_DIR)
  ke=regionSize(THIRD_DIR)
  varCount=regionSize(STRUCTSIZE)
  isFace = (gridDataStruct==FACEX).and.(axis==IAXIS)
  isFace = isFace.or.((gridDataStruct==FACEY).and.(axis==JAXIS))
  isFace = isFace.or.((gridDataStruct==FACEZ).and.(axis==KAXIS))

  
!!  print*,'in applyBcRegion ',varCount,gridDataStruct,guard,axis,face

  do ivar = 1,varCount
     if(mask(ivar)) then
        call gr_bcMapBcType(bcTypeActual,bcType,ivar,gridDataStruct,axis,face,idest)
        sign = 1
        if (gridDataStruct==CENTER) then
#ifdef VELX_VAR
           if ((axis==IAXIS).and.(ivar==VELX_VAR))sign=-1
#endif
#ifdef VELY_VAR
           if((axis==JAXIS).and.(ivar==VELY_VAR))sign=-1
#endif
#ifdef VELZ_VAR
           if((axis==KAXIS).and.(ivar==VELZ_VAR))sign=-1
#endif
        end if
        
        
        if(face==LOW) then
           select case (bcTypeActual)
           case(REFLECTING)
              k = 2*guard+1
              if(isFace)k=k+1
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)*real(sign)
              end do
           case(USER_DEFINED)
              if (((ivar .GE. SPECIES_BEGIN .AND. ivar .LE. SPECIES_END) .OR. &
                   ivar == NI56_MSCALAR .OR. &
                   ivar == NUMP_VAR) .AND. gridDataStruct==CENTER) then
                 select case (ivar)
                 case(C_SPEC,O_SPEC)
                    fixVal = 0.5
                 case(NUMP_VAR)
                    fixVal = 1.e-29
                 case default
                    fixVal = 0.0
                 end select
                 regionData(1:guard,1:je,1:ke,ivar)= fixVal
              else if (( &
                   ivar == VELX_VAR .OR. &
                   ivar == VELY_VAR) .AND. gridDataStruct==CENTER) then
                 select case (ivar)
                 case(VELX_VAR,VELY_VAR)
                    do k=1,ke
                       do j=1,je
                          fixVal = minval(regionData(guard+1:2*guard,j,k,ivar))
                          if (sign < 0.0) then
                             fixVal = minval(regionData(guard+1:2*guard,j,k,ivar))
                          else if (ivar==VELX_VAR .AND. &
                               maxval(regionData(guard+1:2*guard,j,k,ivar)) .GE. 0.0) then
                             fixVal = maxval(regionData(guard+1:2*guard,j,k,ivar))
                          else
                             fixVal = minval(regionData(guard+1:2*guard,j,k,ivar))
                          end if
                          regionData(1:guard,j,k,ivar)= fixVal
                       end do
                    end do
                 end select
              else              !ivar
                 k = 2*guard+1
                 if(isFace)k=k+1
                 do i = 1,guard
                    regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)
                 end do
              end if

           case(DIRICHLET)
              do i = 1,guard
                 regionData(guard+1-i,1:je,1:ke,ivar)= (1-2*i)*regionData(guard+1,1:je,1:ke,ivar)
              end do
           case(GRIDBC_MG_EXTRAPOLATE)
              do i = 1,guard
                 regionData(guard+1-i,1:je,1:ke,ivar)= (1+i)*regionData(guard+1,1:je,1:ke,ivar) &
                      - i*regionData(guard+2,1:je,1:ke,ivar)
              end do

           case(OUTFLOW,HYDROSTATIC_F2_NVOUT,HYDROSTATIC_F2_NVDIODE,HYDROSTATIC_F2_NVREFL, &
                        HYDROSTATIC_NVOUT,HYDROSTATIC_NVDIODE,HYDROSTATIC_NVREFL)
!!              print*,'since face was low',je,ke,ivar
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(guard+1,1:je,1:ke,ivar)
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
           case default
!!              print*,'boundary is',bcType
!!              call Driver_abortFlash("unsupported boundary condition on Lower Face")
           end select
           
        else  !(face==HIGH)
           
           select case (bcTypeActual)
           case(REFLECTING)
              k = 2*guard+1
              if(isFace)k=k+1
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)*real(sign)
              end do
           case(USER_DEFINED)
              if (((ivar .GE. SPECIES_BEGIN .AND. ivar .LE. SPECIES_END) .OR. &
                   ivar == NI56_MSCALAR .OR. &
                   ivar == NUMP_VAR) .AND. gridDataStruct==CENTER) then
                 select case (ivar)
                 case(C_SPEC,O_SPEC)
                    fixVal = 0.5
                 case(NUMP_VAR)
                    fixVal = 1.e-29
                 case default
                    fixVal = 0.0
                 end select
                 regionData(guard+1:2*guard,1:je,1:ke,ivar)= fixVal
              else if (( &
                   ivar == VELX_VAR .OR. &
                   ivar == VELY_VAR) .AND. gridDataStruct==CENTER) then
                 select case (ivar)
                 case(VELX_VAR,VELY_VAR)
                    do k=1,ke
                       do j=1,je
                          if (sign < 0.0) then
                             fixVal = maxval(regionData(1:guard,j,k,ivar))
                          else if (ivar==VELY_VAR .AND. &
                               minval(regionData(1:guard,j,k,ivar)) < 0.0) then
                             fixVal = minval(regionData(1:guard,j,k,ivar))
                          else
                             fixVal = maxval(regionData(1:guard,j,k,ivar))
                          end if
                          regionData(guard+1:2*guard,j,k,ivar)= fixVal
                       end do
                    end do
                 end select
              else              !ivar
                 k = 2*guard+1
                 if(isFace)k=k+1
                 do i = 1,guard
                    regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)
                 end do
              end if

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

           case(OUTFLOW,HYDROSTATIC_F2_NVOUT,HYDROSTATIC_F2_NVDIODE,HYDROSTATIC_F2_NVREFL, &
                        HYDROSTATIC_NVOUT,HYDROSTATIC_NVDIODE,HYDROSTATIC_NVREFL)
              k=guard
              if(isFace)k=k+1
!!              print*,'since face was high',k,je,ke,ivar
              do i = 1,guard
                 regionData(k+i,1:je,1:ke,ivar)= regionData(k,1:je,1:ke,ivar)
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
           case default
!!              print*,'boundary is',bcType
!!              call Driver_abortFlash("unsupported boundary condition on Upper Face")
           end select
        end if
     end if
  end do

  return
end subroutine Grid_bcApplyToRegion
