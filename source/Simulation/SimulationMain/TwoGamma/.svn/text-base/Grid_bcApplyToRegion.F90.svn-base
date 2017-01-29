!!****if* source/Simulation/SimulationMain/TwoGamma/Grid_bcApplyToRegion
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
!!                                       integer(IN)  :: blockHandle,
!!                                       integer(IN)  :: secondDir,
!!                                       integer(IN)  :: thirdDir,
!!                                       integer(IN)  :: endPoints(LOW:HIGH,MDIM),
!!                                       integer(IN)  :: blkLimitsGC(LOW:HIGH,MDIM),
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
!!  This routine applies the boundary conditions on a given face (lowerface
!!  or upperface) along a given axis, by using and setting values
!!  for all variables in the gridDataStruct that are not masked out. The 
!!  argument "mask" has the information about the masked variables.
!! 
!!   If (face=LOW)
!!     For all iVar where mask(iVar) is TRUE:
!!       regionData(1:guard,:,:,iVar) =  boundary values
!!   If (face=HIGH)
!!     For all iVar where mask(iVar) is TRUE:
!!       regionData(regionSize(BC_DIR)-guard+1:regionSize(BC_DIR),:,:,iVar) =  boundary values
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
!!
!!  This routine supports only simple boundary conditions that are applied strictly
!!  directionally, and have no need for other grid information such as the coordinates etc.
!!  Currently supported boundary conditions are "OUTFLOW", "REFLECTING" and "DIODE".
!!  The "PERIODIC" boundary conditions are automatically applied in the process of filling
!!  the guard cells all over the domain, and therefore do not need to call this routine
!!
!!  This implementation specific to the TwoGamma simulation additionally supports
!!  a USER_DEFINED boundary condition for matter flowing in from the left (or lower) face.
!!
!!  If the user wishes to apply different boundary conditions, they can either
!!  use the interface Grid_bcApplyToRegionSpecialized, or make a copy of this
!!  routine in their simulation directory and customize it.
!!
!!  
!! ARGUMENTS 
!!
!!  bcType - the type of boundary condition being applied
!!  gridDataStruct - the Grid dataStructure
!!  guard -    number of guardcells 
!!  axis  - the dimension along which to apply boundary conditions,
!!          can take values of IAXIS, JAXIS and KAXIS
!!  face    -  can take values LOW and HIGH, defined in constants.h,
!!             to indicate whether to apply boundary on lowerface or 
!!             upperface
!!  regionData     : the extracted region from a block of permanent storage of the 
!!                   specified data structure. Its size is given by regionSize.
!!  regionSize     : regionSize(BC_DIR) contains the size of each row of data
!!                   in the regionData array.  For the common case of guard=4,
!!                   regionSize(BC_DIR) will be 8 for cell-centered data structures
!!                   (e.g., when gridDataStruct=CENTER) and either 8 or 9 for face-
!!                   centered data, depending on the direction given by axis.
!!                   regionSize(SECOND_DIR) contains the number of rows along the second
!!                   direction, and regionSize(THIRD_DIR) has the number of rows
!!                   along the third direction. regionSize(GRID_DATASTRUCT) contains the
!!                   number of variables in the data structure.
!!  mask - if present applies boundary conditions to only selected variables
!!  applied - is set true if this routine has handled the given bcType, otherwise it is 
!!            set to false.
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
!!  idest   - ignored in this implementation specific to the TwoGamma setup.
!!
!!  NOTES 
!!
!!            This routine can be used with all the mesh packages supported.
!!            The mesh packages extract the small vectors relevant to
!!            boundary conditions calculations from their Grid data 
!!            structures. If users wish to apply a different boundary condition, 
!!            they should look at routine Grid_bcApplyToRegionSpecialized.
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

  use Simulation_data, ONLY : sim_rho1,sim_rho2,sim_temp1,sim_temp2,&
                              sim_gammac1,sim_gammac2,&
                              sim_p0,sim_cvelx,sim_int1,sim_small
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY: gr_meshMe   ! only needed for print* in error case

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

  logical,dimension(regionSize(STRUCTSIZE)) :: localMask
  integer :: i,j, k,ivar,ind,sign,je,ke,n,varCount
  logical :: isFace

  if (gridDataStruct/=CENTER) then
     print*,'boundary is',bcType,gr_meshMe,face,gridDataStruct
     call Driver_abortFlash("[Grid_applyBCEdge] Simulation TwoGamma does not support face variables")
  end if

  je=regionSize(SECOND_DIR)
  ke=regionSize(THIRD_DIR)
  varCount=regionSize(STRUCTSIZE)
  isFace = (gridDataStruct==FACEX).and.(axis==IAXIS)
  isFace = isFace.or.((gridDataStruct==FACEY).and.(axis==JAXIS))
  isFace = isFace.or.((gridDataStruct==FACEZ).and.(axis==KAXIS))

!!  print*,'in Grid_bcApplyToRegion ',varCount,gridDataStruct,guard,axis,face
  do ivar = 1,varCount
     if(mask(ivar)) then
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
           select case (bcType)
           case(REFLECTING)
              k = 2*guard+1
              if(isFace)k=k+1
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)*sign
              end do
        
           case(OUTFLOW)
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
           case(USER_DEFINED)
              select case(ivar)
              case(GAMC_VAR)
                 regionData(1:guard,1:je,1:ke,ivar)=sim_gammac1
              case(DENS_VAR)
                 regionData(1:guard,1:je,1:ke,ivar)=sim_rho1
              case(PRES_VAR)
                 regionData(1:guard,1:je,1:ke,ivar)=sim_p0
              case(VELX_VAR)
                 regionData(1:guard,1:je,1:ke,ivar)=sim_cvelx
              case(VELY_VAR)
                 regionData(1:guard,1:je,1:ke,ivar)=0.0
              case(VELZ_VAR)
                 regionData(1:guard,1:je,1:ke,ivar)=0.0
              case(ENER_VAR)
                 regionData(1:guard,1:je,1:ke,ivar)=max(0.5*(sim_cvelx**2)+sim_int1,sim_small)
              case(EINT_VAR)
                 regionData(1:guard,1:je,1:ke,ivar)=max(sim_int1,sim_small)
              case(SPECIES_BEGIN)
                 regionData(1:guard,1:je,1:ke,ivar)=1.0e0-sim_small
              case(SPECIES_BEGIN+1)
                 regionData(1:guard,1:je,1:ke,ivar)=sim_small
              case(TEMP_VAR)
                 regionData(1:guard,1:je,1:ke,ivar)=sim_temp1
#ifdef TION_VAR
              case(TION_VAR,TELE_VAR)
                 regionData(1:guard,1:je,1:ke,ivar)=sim_temp1
#endif
#ifdef ERAD_VAR
              case(ERAD_VAR)
                 regionData(1:guard,1:je,1:ke,ivar)=0.0
#endif
#ifdef PRAD_VAR
              case(PRAD_VAR)
                 regionData(1:guard,1:je,1:ke,ivar)=0.0
#endif
#ifdef TRAD_VAR
              case(TRAD_VAR)
                 regionData(1:guard,1:je,1:ke,ivar)=0.0
#endif
              case default      !like OUTFLOW
                 do i = 1,guard
                    regionData(i,1:je,1:ke,ivar)= regionData(guard+1,1:je,1:ke,ivar)
                 end do
              end select
           case default
              print*,'boundary is',bcType
              call Driver_abortFlash("unsupported boundary condition on Lower Face")
           end select
           
        else
           
           select case (bcType)
           case(REFLECTING)
              k = 2*guard+1
              if(isFace)k=k+1
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)*sign
              end do
              
           case(OUTFLOW)
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
           case(USER_DEFINED)
              print*,'boundary is',bcType,gr_meshMe,face
              call Driver_abortFlash("Simulation TwoGamma does not support USER_DEFINED boundary on Upper Face")
           case default
              print*,'boundary is',bcType
              call Driver_abortFlash("unsupported boundary condition on Upper Face")
           end select
        end if
     end if
  end do

  applied = .TRUE.              ! if we got to here... - KW
  return
end subroutine Grid_bcApplyToRegion
