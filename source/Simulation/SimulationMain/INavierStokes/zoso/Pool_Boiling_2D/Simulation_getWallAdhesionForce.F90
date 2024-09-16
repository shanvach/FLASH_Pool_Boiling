subroutine Simulation_getWallAdhesionForce(fx_p,fy_p,fx_visc,fy_visc,fx_int,fy_int,fy_g,fx_tr,fy_tr,dt)

      use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_getBlkCenterCoords, Grid_getDeltas
      use IncompNS_data, ONLY : ins_gravX,ins_gravY
      use Driver_data, ONLY: dr_nstep
      use Multiphase_data, ONLY : mph_rho1,mph_rho2,mph_vis1,mph_vis2
                              

#include "Flash.h"
#include "constants.h"


      implicit none
     ! real, dimension(:,:,:),intent(in)    :: u,v,nx,ny,curv
      integer :: blockID
     ! real, intent(in) :: dy,dx
      !integer :: jy1,jy2,ix1,ix2,kz1,kz2
      integer,  dimension(MAXBLOCKS) :: blockList
      integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
      integer :: i,j,k,lb,blockCount
      real del(MDIM),bsize(MDIM),coord(MDIM)
      real, dimension(2,MDIM) :: boundBox
      real :: xcell, zcell,ycell,xcellp
      real :: tol=1E-13,xe = 2.0,ye = 2.0
      
      real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData

      real,    INTENT(IN) :: dt
      real, intent(out) :: fx_p,fy_p,fx_visc,fy_visc,fx_int,fy_int
      real, intent(out) :: fy_g,fx_tr,fy_tr
      real :: dudx,dudy,dvdx,dvdy,dA,dens,u_up,u_lo,v_ri,v_le,un,vn

      !call Grid_getDeltas(blockID,del)
      !call Grid_getBlkCenterCoords(blockId,coord)
      !!call Grid_getBlkBoundBox(blockId,boundBox)

      do lb = 1,blockCount

     blockID = blockList(lb)

     call Grid_getBlkBoundBox(blockId,boundBox)

     bsize(:) = boundBox(2,:) - boundBox(1,:)

     call Grid_getBlkCenterCoords(blockId,coord)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)


     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
      do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           
             xcell  = coord(IAXIS) - bsize(IAXIS)/2.0 +  &
                    real(i - NGUARD - 1)*del(IAXIS)  +  &
                    0.5*del(IAXIS)
             ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                    real(j - NGUARD - 1)*del(JAXIS)  +  &
                    0.5*del(JAXIS)
             xcellp  = coord(IAXIS) - bsize(IAXIS)/2.0 +  &
                    real(i+1 - NGUARD - 1)*del(IAXIS)  +  &
                    0.5*del(IAXIS)

            if (solnData(DFUN_VAR,i,j,1) .gt. 0 .and. solnData(DFUN_VAR,i+1,j,1) .lt. 0) then
	     	dA = abs(solnData(NRMX_VAR,i+1,j,1))*del(JAXIS) + abs(solnData(NRMY_VAR,i+1,j,1))*del(IAXIS)
		fx_p = fx_p +  solnData(PRES_VAR,i+1,j,1)*solnData(NRMX_VAR,i+1,j,1)*dA
		fy_p = fy_p +  solnData(PRES_VAR,i+1,j,1)*solnData(NRMY_VAR,i+1,j,1)*dA	
	    

                un =  (facexData(VELC_FACE_VAR,i,j,1)+facexData(VELC_FACE_VAR,i+1,j,1))*0.5
                vn =  (faceyData(VELC_FACE_VAR,i,j,1)+faceyData(VELC_FACE_VAR,i-1,j,1)+faceyData(VELC_FACE_VAR,i,j,1)+faceyData(VELC_FACE_VAR,i,j-1,1))*0.25

	        fx_int = fx_int + (0.5*(mph_rho1 + mph_rho2))*(un*-solnData(NRMX_VAR,i+1,j,1) + vn*-solnData(NRMY_VAR,i+1,j,1))*un*dA
		fy_int = fy_int + (0.5*(mph_rho1 + mph_rho2))*(un*-solnData(NRMX_VAR,i+1,j,1) + vn*-solnData(NRMY_VAR,i+1,j,1))*vn*dA
	

                dudx =   (facexData(VELC_FACE_VAR,i+1,j,1)-facexData(VELC_FACE_VAR,i,j,1))/del(IAXIS)
                u_up =  (facexData(VELC_FACE_VAR,i+1,j+1,1)+facexData(VELC_FACE_VAR,i+1,j,1)+facexData(VELC_FACE_VAR,i,j,1)+facexData(VELC_FACE_VAR,i,j+1,1))*0.25
	        u_lo =  (facexData(VELC_FACE_VAR,i+1,j-1,1)+facexData(VELC_FACE_VAR,i+1,j,1)+facexData(VELC_FACE_VAR,i,j,1)+facexData(VELC_FACE_VAR,i,j-1,1))*0.25
		dudy = (u_up - u_lo)/del(JAXIS)
                v_ri =  (faceyData(VELC_FACE_VAR,i+1,j-1,1)+faceyData(VELC_FACE_VAR,i+1,j,1)+faceyData(VELC_FACE_VAR,i,j,1)+faceyData(VELC_FACE_VAR,i,j-1,1))*0.25
                v_le =  (faceyData(VELC_FACE_VAR,i-1,j-1,1)+faceyData(VELC_FACE_VAR,i-1,j,1)+faceyData(VELC_FACE_VAR,i,j,1)+faceyData(VELC_FACE_VAR,i,j-1,1))*0.25
                dvdx = (v_ri - v_le)/del(IAXIS)
                dvdy =  (faceyData(VELC_FACE_VAR,i+1,j,1)-faceyData(VELC_FACE_VAR,i+1,j-1,1))/del(JAXIS)
               
                fx_visc = fx_visc + (mph_vis2)*(dudx*-solnData(NRMX_VAR,i+1,j,1) + 0.5*(dudy+dvdx)*-solnData(NRMY_VAR,i+1,j,1))*dA
                fy_visc = fy_visc + (mph_vis2)*(dvdy*-solnData(NRMY_VAR,i+1,j,1) + 0.5*(dudy+dvdx)*-solnData(NRMX_VAR,i+1,j,1))*dA
            
	   end if
             
               fy_g = fy_g+ solnData(PFUN_VAR,i,j,1)*(del(IAXIS))*(del(JAXIS))*(mph_rho1)*(ins_gravY)*9.81
	       fx_tr = fx_tr + (mph_rho1)*solnData(PFUN_VAR,i,j,1)*(del(IAXIS))*(del(JAXIS))*(0.5*(facexData(VELC_FACE_VAR,i-1,j,1)+facexData(VELC_FACE_VAR,i,j,1)))
               fy_tr = fy_tr + (mph_rho1)*solnData(PFUN_VAR,i,j,1)*(del(IAXIS))*(del(JAXIS))*(0.5*(faceyData(VELC_FACE_VAR,i,j-1,1)+faceyData(VELC_FACE_VAR,i,j,1)))
    


     end do
   end do


end do
end subroutine
