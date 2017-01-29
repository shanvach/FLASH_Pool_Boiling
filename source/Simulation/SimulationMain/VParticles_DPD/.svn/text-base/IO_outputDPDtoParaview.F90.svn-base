subroutine IO_outputDPDtoParaview(filename,mype,time,dt,istep,count, &
     blockList,blockCount,firstfileflag,particles,p_count)
  
#include "constants.h"
#include "Flash.h"

  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_getBlkIndexLimits, &
       Grid_getBlkBoundBox,Grid_getBlkCenterCoords
  
  use HDF5  
#ifdef FLASH_GRID_PARAMESH
  use physicaldata, only : interp_mask_unk, interp_mask_unk_res
#endif


  implicit none
#include "Flash_mpi.h"
#include "Particles.h"
#include "io_flash.h" 
 
  integer, intent(in)   :: mype,istep,count,firstfileflag
  integer, intent(in)   :: blockCount
  integer, intent(in)   :: blockList(MAXBLOCKS)
  real, intent(in)      :: time,dt
  character(len=100),intent(in)  :: filename
  integer,INTENT(IN):: p_count
  real,dimension(NPART_PROPS,p_count),INTENT(IN)::particles
  
  ! Local variables    
  integer :: numblocks,var,i,j,k,lb,nxc,nyc,nzc
  character(6) :: index_lb,index_mype

  real:: xedge(NXB+1),xcell(NXB+1)
  real:: yedge(NYB+1),ycell(NYB+1)
  real:: zedge(NZB+1),zcell(NZB+1)
  real:: intsx(NXB+1), intsy(NYB+1), intsz(NZB+1)

  integer ::  blockID,bodycounter

  real ::  del(3),dx,dy,dz
  real, dimension(MDIM)  :: coord,bsize
  real ::  boundBox(2,MDIM)

  ! HDF5 variables
  integer :: h5err, h5_outfile_id,group_id
  integer(HID_T) :: dset_id, space_id
  integer(HSIZE_T), DIMENSION(3) :: dims ! size read/write buffer


  character(len=100) :: filenameh5,bodyname,gridName,filenamexmf
 

  ! ---
  ! --- Initialize the Fortran interface 
  ! --- 
  CALL h5open_f(h5err)


  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1
  nzc = NZB + NGUARD + 1

  intsx    = (/ (real(i), i=0,NXB) /)
  intsy    = (/ (real(i), i=0,NYB) /)
  intsz    = (/ (real(i), i=0,NZB) /)

  do lb = 1,blockcount
     
     blockID =  blockList(lb)
     
     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)
     dx = del(IAXIS)
     dy = del(JAXIS)
     dz = del(KAXIS)
     
     ! Get Coord and Bsize for the block:
     ! Bounding box:
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)

     call Grid_getBlkCenterCoords(blockId,coord)
     
     xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
     xcell = xedge(:) + dx/2.0;
     
     yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
     ycell = yedge(:) + dy/2.0;
     
     zedge = coord(KAXIS) - bsize(KAXIS)/2.0 + dz*intsz;
     zcell = zedge(:) + dz/2.0;
         
     
    if (lb.eq.1) then
        ! HDF5
        ! write solution data to DPD_parts_.XXXX.XXXX
       write(*,*) '  '
       write(filenameh5,'(A14,I4.4,".",I4.4,".h5")') TRIM(filename),count,mype
       write(*,*) 'The hdf5 file is ',filenameh5
       call H5Fcreate_f(trim(adjustl(filenameh5)),H5F_ACC_TRUNC_F, h5_outfile_id, h5err)
       call dpd_io_checkHDFerr(h5err, 'FAILURE TO OPEN FILE')
       group_id = h5_outfile_id
     
     
       ! XDMF
       write(filenamexmf,'(A,I4.4,".",I4.4,".xmf")') TRIM(filename),count,mype
       open(unit=30,file=filenamexmf,STATUS='UNKNOWN')
       write(30,'(A)') '<?xml version="1.0" ?>'
       write(30,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
       write(30,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2">'
       write(30,'(A,I7,A)')' <Information Name="Iteration" Value="',istep,'"/>'
       write(30,'(A,f6.4,A)') ' <Time TimeType="Single" Value="',time,'"/>'
       write(30,'(A)') '  <Domain>'
       write(30,'(A)') '   <Grid GridType="Collection" CollectionType="Spatial">'
    end if
    

     write(30,'(A,I3,A)') '    <Grid Name="Fluidmesh[',lb-1,']">'
     write(30,'(A,I4,I4,I4,A)') '     <Topology TopologyType="3DRectMesh" NumberOfElements="',NZB+1,NYB+1,NXB+1,'"/>'
     write(30,'(A)') '      <Geometry GeometryType="VXVYVZ">'
     
     ! Write the x
     dims=(/(NXB+1),1,1/)
     
     write(bodyname,'("X_",I3.3)')lb
     call H5Screate_simple_f(1, dims, space_id, h5err)
     call H5Dcreate_f(group_id,TRIM(bodyname), H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
     call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, xedge, dims, h5err)
     call H5Dclose_f(dset_id, h5err)
     call H5Sclose_f(space_id, h5err)
     write(30,'(A,I4,A)') '       <DataItem Dimensions="',(NXB+1),'" NumberType="Float" Precision="8" Format="HDF">'
     write(30,'(A,A,A,A,A)') '        ',TRIM(filenameh5),':/',TRIM(bodyname)
     write(30,'(A)') '       </DataItem>'
     
     write(bodyname,'("Y_",I3.3)')lb
     dims=(/(NYB+1),1,1/)
     call H5Screate_simple_f(1, dims, space_id, h5err)
     call H5Dcreate_f(group_id,TRIM(bodyname), H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
     call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, yedge, dims, h5err)
     call H5Dclose_f(dset_id, h5err)
     write(30,'(A,I4,A)') '       <DataItem Dimensions="',(NYB+1),'" NumberType="Float" Precision="8" Format="HDF">'
     write(30,'(A,A,A,A,A)')'         ', TRIM(filenameh5),':/',TRIM(bodyname)
     write(30,'(A)') '       </DataItem>'
     ! write the z
     dims=(/(NZB+1),1,1/)
     write(bodyname,'("Z_",I3.3)')lb
     call H5Screate_simple_f(1, dims, space_id, h5err)
     call H5Dcreate_f(group_id,TRIM(bodyname), H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
     call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, zedge, dims, h5err)
     call H5Dclose_f(dset_id, h5err)
     write(30,'(A,I4,A)') '       <DataItem Dimensions="',NZB+1,'" NumberType="Float" Precision="8" Format="HDF">'
     write(30,'(A,A,A,A,A)')'         ',TRIM(filenameh5),':/',TRIM(bodyname)
     write(30,'(A)') '       </DataItem>'
     write(30,'(A)')'      </Geometry>'
     write(30,'(A)') '    </Grid>'
  end do
     
   write(bodyname,'("parts_",I3.3)')mype
   dims=(/NDIM,p_count,1/)
   call H5Screate_simple_f(2, dims, space_id, h5err)
   call H5Dcreate_f(group_id,TRIM(bodyname), H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
   call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, particles(POSX_PART_PROP:POSZ_PART_PROP,1:p_count), dims, h5err)
   call H5Dclose_f(dset_id, h5err)
   
   write(30,'(A,I3.3,A)') '   <Grid Name="dpd_particles_',mype,'" GridType="Uniform">'
   write(30,'(A,F6.2,A)') '    <Time Value="',time,'" />'
   write(30,'(A,I6,A)')'     <Topology TopologyType="Polyvertex" NodesPerElement="',p_count,'">'
   write(30,'(A)')'      </Topology>'
   write(30,'(A)')'       <Geometry GeometryType="XYZ">'
   write(30,'(A,I6,A)')'        <DataItem DataType="Float" Dimensions="',p_count ,' 3" Format="HDF">'
   write(30,'(A,A,A,A)')'         ', TRIM(filenameh5),':/',TRIM(bodyname)
   write(30,'(A)')'       </DataItem>'
   write(30,'(A)')'       </Geometry>'
   write(30,'(A)')'    </Grid>'
 
   !write(*,*) 'blockcount',blockcount
   ! Close the group
   if (blockcount.ge.1) then
      write(30,'(A)') '   </Grid>'
      write(30,'(A)') '  </Domain>'
      write(30,'(A)') '</Xdmf>'
      close(30)
      if( group_id /= h5_outfile_id ) then
         call H5Gclose_f (group_id, h5err)
      end if
   end if
   
   !
   ! close the file
   !
   call H5Fclose_f(h5_outfile_id, h5err)
   call dpd_io_checkHDFerr(h5err,'failure to close HDF5 file')
   ! --- 
   ! Clean up- HDF5 stuf
   CALL H5close_f(h5err)
   call dpd_io_checkHDFerr(h5err,'failure to clean-up HDF5 stuff')
 end subroutine IO_outputDPDtoParaview
 
