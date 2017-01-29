interface Grid_getCoordPtr

   subroutine Grid_getCoordPtr(blockId, axis, coordPtr)
     integer, intent(in) :: blockId, axis
     real, dimension(:,:), pointer :: coordPtr
   end subroutine Grid_getCoordPtr
end interface

interface Grid_getBlkPtr
subroutine Grid_getBlkPtr(blockId, dataPtr)
integer, intent(in) :: blockId
real,dimension(:,:,:,:), pointer :: dataPtr
end subroutine Grid_getBlkPtr
end interface

interface Grid_releaseBlockPtr
  subroutine Grid_releaseBlockPtr(dataPtr)
     real, pointer :: dataPtr(:,:,:,:)
  end subroutine Grid_releaseBlockPtr
end interface

interface Grid_releaseCoordPtr
  subroutine Grid_releaseCoordPtr(dataPtr)
     real, pointer :: dataPtr(:,:)
  end subroutine Grid_releaseCoordPtr
end interface

interface Grid_getSavedVarPtr
subroutine Grid_getSavedVarPtr(blockId, dataPtr)
implicit none
integer, intent(in) :: blockID
real,dimension(:,:,:,:), pointer :: dataPtr
end subroutine Grid_getSavedVarPtr
end interface

interface Grid_releaseSavedVarPtr
  subroutine Grid_releaseSavedVarPtr(dataPtr)
     real, pointer :: dataPtr(:,:,:,:)
  end subroutine Grid_releaseSavedVarPtr
end interface

