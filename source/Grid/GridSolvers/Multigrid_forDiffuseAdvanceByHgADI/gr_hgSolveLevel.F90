!!****if* source/Grid/GridSolvers/Multigrid_forDiffuseAdvanceByHgADI/gr_hgSolveLevel
!!
!! NAME
!!  gr_hgSolveLevel
!!
!! SYNOPSIS
!!  call gr_hgSolveLevel(integer, intent(in) :: level,
!!                      integer, intent(in) :: gr_iSource,
!!                      integer, intent(in) :: gr_iSoln,
!!                      external            :: SolveBlock,
!!                      integer, intent(in) :: LeafFlag,
!!                      OPTIONAL,real(IN)   :: dt,
!!                      OPTIONAL,real(IN)   :: chi,
!!                      OPTIONAL,real(IN)   :: theta)
!! 
!! DESCRIPTION
!!
!!  Block-per-block solve Diffusion on all the blocks at level using boundary
!!  conditions on the exterior faces.
!! 
!! ARGUMENTS
!!  
!!  level        - the level to solve on
!!  gr_iSource   - the grid variable holding the source term
!!  gr_iSoln     - the grid variable in which the solution resides
!!  SolveBlock   - a function that can solve (or smooth) locally
!!  LeafFlag     - 0 => all blocks on the level solved
!!                 1 => only leaf blocks solved
!!                 2 => only parent blocks solved
!!  dt         - time step, to be passed down 
!!  chi        - a factor, to be passed down 
!!  theta      - 0.5 => Crank Nicholson Implicit Time Integration.
!!               0.0 => Explicit Time Integration.
!!
!! NOTES
!!
!!  We don't have anything equivalent to the Jacobi sweeps that are used
!!  when solving the Poisson problem. Jacobi sweeps code is removed.
!!
!!***

!!REORDER(5): unk

subroutine gr_hgSolveLevel(level, gr_iSource, gr_iSoln, SolveBlock, LeafFlag, dt, chi,theta)

!==================================================================

  use gr_hgData, ONLY: gr_hgBndTypes, &
       hg_ili, hg_iui, hg_jli, hg_jui, hg_kli, hg_kui
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getDeltas
  use workspace, ONLY : work
  use tree, ONLY : lnblocks,lrefine,nodetype,bsize
  use physicaldata, ONLY : unk
  
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use gr_hgInterface, ONLY: gr_hgBndry

  implicit none
#include "Flash.h"
#include "constants.h"
#include "Multigrid.h"
#include "Flash_mpi.h"

  integer, intent(in) :: gr_iSource, gr_iSoln, level, LeafFlag
  real,intent(IN),OPTIONAL :: dt, chi, theta

  external               SolveBlock
  
  integer                      :: b, i, j, k, n
  integer                      :: nblockssolved

  integer,parameter            :: xtraSolnCellsI=1,xtraSolnCellsJ=(1*K2D),xtraSolnCellsK=(1*K3D)
  integer,parameter            :: hgSL_ile = 1-xtraSolnCellsI, hgSL_iue = NXB+xtraSolnCellsI
  integer,parameter            :: hgSL_jle = 1-xtraSolnCellsJ, hgSL_jue = NYB+xtraSolnCellsJ
  integer,parameter            :: hgSL_kle = 1-xtraSolnCellsK, hgSL_kue = NZB+xtraSolnCellsK
  integer,parameter            :: xtraSrcCellsI=0, xtraSrcCellsJ=(0*K2D), xtraSrcCellsK=(0*K3D)
  integer,parameter            :: hgSR_ile = 1-xtraSrcCellsI, hgSR_iue = NXB+xtraSrcCellsI
  integer,parameter            :: hgSR_jle = 1-xtraSrcCellsJ, hgSR_jue = NYB+xtraSrcCellsJ
  integer,parameter            :: hgSR_kle = 1-xtraSrcCellsK, hgSR_kue = NZB+xtraSrcCellsK

  real :: soln(1-xtraSolnCellsI:NXB+xtraSolnCellsI,1-xtraSolnCellsJ:NYB+xtraSolnCellsJ,1-xtraSolnCellsK:NZB+xtraSolnCellsK)
  real :: oldSrc(1-xtraSrcCellsI:NXB+xtraSrcCellsI,1-xtraSrcCellsJ:NYB+xtraSrcCellsJ,1-xtraSrcCellsK:NZB+xtraSrcCellsK)

  real :: rhs(NXB,NYB,NZB)

  real, dimension(MDIM)        :: deltas

  logical                      :: SolveThisBlock
  integer                      :: bnd_type

  character                    :: level_timer_label*40

  !====================================================================


  write(unit=level_timer_label, FMT='(A11,I2)') "solvelevel", level
  if (not(present(dt))) call Driver_abortFlash('This gr_hgSolveLevel implementation requires dt to be present')
  if (not(present(chi))) call Driver_abortFlash('This gr_hgSolveLevel implementation requires chi to be present')
  call Timers_start(level_timer_label)

  nblockssolved = 0

  do b = 1, lnblocks
     
     SolveThisBlock = (lrefine(b) == level)
     if (LeafFlag == MG_NODES_LEAF_ONLY) then
        SolveThisBlock = (SolveThisBlock .and. (nodetype(b) == LEAF))
     else if (LeafFlag == MG_NODES_PARENT_ONLY) then
        SolveThisBlock = (SolveThisBlock .and. (nodetype(b) == PARENT_BLK))
     endif

     if (SolveThisBlock) then
        nblockssolved = nblockssolved+1

        !Let's try a standard way of finding delta x rather than this convolution        
        call Grid_getDeltas(b,deltas) 

        do k = 1, NZB
           do j = 1, NYB
              do i = 1, NXB
                 soln(i,j,k) = unk(gr_iSoln,i+NGUARD,j+K2D*NGUARD,k+K3D*NGUARD,b)                 
              enddo
           enddo
        enddo
         do k = 1, NZB
           do j = 1, NYB
              do i = 1-xtraSrcCellsI, NXB+xtraSrcCellsI
                 oldSrc(i,j,k) = unk(gr_iSource,i+NGUARD,j+K2D*NGUARD,k+K3D*NGUARD,b)
              enddo
           enddo
        enddo


        bnd_type = gr_hgBndTypes(2*NDIM-1) !low face highest dimension, should really check all directions...
        if (bnd_type == MG_BND_PERIODIC) then
           ! Only the coarsest mesh level is treated as periodic; the rest obtain
           ! boundary values by interpolation and are treated as Dirichlet.
           if (level == 1) then
              bnd_type = MG_BND_PERIODIC
           else
              bnd_type = MG_BND_DIRICHLET
           endif
        else if ((bnd_type == MG_BND_DIRICHLET) .or. &
                 (bnd_type == MG_BND_GIVENVAL)) then
           bnd_type = MG_BND_DIRICHLET
        else if ((bnd_type == MG_BND_NEUMANN)) then
           bnd_type = MG_BND_DIRICHLET
        else
           call Driver_abortFlash("gr_hgSolveLevel found an unrecognized bnd_type!")
        endif
        
        ! Update the solution function GC to account for given-value boundary conditions.
        if (bnd_type == MG_BND_PERIODIC) then
           

           ! X - Direction BC.
           soln(0,    1:NYB,1:NZB) = unk(gr_iSoln,NGUARD+NXB+1, 1+NGUARD*K2D:NYB+NGUARD*K2D, 1+NGUARD*K3D:NZB+NGUARD*K3D,b)
           soln(NXB+1,1:NYB,1:NZB) = unk(gr_iSoln,NGUARD      , 1+NGUARD*K2D:NYB+NGUARD*K2D, 1+NGUARD*K3D:NZB+NGUARD*K3D,b)

           ! Y - Direction BC.
           if (NDIM .ge. 2) then
               soln(1:NXB,0    ,1:NZB) = unk(gr_iSoln,1+NGUARD:NXB+NGUARD, NGUARD*K2D+NYB+1      , 1+NGUARD*K3D:NZB+NGUARD*K3D,b)
               soln(1:NXB,NYB+1,1:NZB) = unk(gr_iSoln,1+NGUARD:NXB+NGUARD, NGUARD*K2D, 1+NGUARD*K3D:NZB+NGUARD*K3D,b)
           endif

           ! Z - Direction BC.
           if (NDIM .eq. 3) then
               soln(1:NXB,1:NYB,0)     = unk(gr_iSoln,1+NGUARD:NXB+NGUARD, 1+NGUARD*K2D:NYB+NGUARD*K2D, NGUARD*K3D,b)
               soln(1:NXB,1:NYB,NZB+1) = unk(gr_iSoln,1+NGUARD:NXB+NGUARD, 1+NGUARD*K2D:NYB+NGUARD*K2D, 1+NZB+NGUARD*K3D,b)
           endif
           
          
        endif 

        if (bnd_type == MG_BND_DIRICHLET) then

           ! X - Direction BC.
           soln(0,    1:NYB,1:NZB) = unk(gr_iSoln,NGUARD      , 1+NGUARD*K2D:NYB+NGUARD*K2D, 1+NGUARD*K3D:NZB+NGUARD*K3D,b)
           soln(NXB+1,1:NYB,1:NZB) = unk(gr_iSoln,NGUARD+NXB+1, 1+NGUARD*K2D:NYB+NGUARD*K2D, 1+NGUARD*K3D:NZB+NGUARD*K3D,b)

           ! Y - Direction BC.
           if (NDIM .ge. 2) then
               soln(1:NXB,0    ,1:NZB) = unk(gr_iSoln,1+NGUARD:NXB+NGUARD, NGUARD*K2D      , 1+NGUARD*K3D:NZB+NGUARD*K3D,b)
               soln(1:NXB,NYB+1,1:NZB) = unk(gr_iSoln,1+NGUARD:NXB+NGUARD, 1+NYB+NGUARD*K2D, 1+NGUARD*K3D:NZB+NGUARD*K3D,b)
           end if

           ! Z - Direction BC
           if (NDIM .eq. 3) then
               soln(1:NXB,1:NYB,0)     = unk(gr_iSoln,1+NGUARD:NXB+NGUARD, 1+NGUARD*K2D:NYB+NGUARD*K2D, NGUARD*K3D,b)
               soln(1:NXB,1:NYB,NZB+1) = unk(gr_iSoln,1+NGUARD:NXB+NGUARD, 1+NGUARD*K2D:NYB+NGUARD*K2D, 1+NZB+NGUARD*K3D,b)
           endif

        endif

        call Timers_start("SolveBlock")
        !
        call SolveBlock (soln, NXB, NYB, NZB, &
                         deltas(IAXIS), deltas(JAXIS), deltas(KAXIS), bnd_type, level, &
                         dt, chi, theta, oldSrc,&
                         hgSL_ile,hgSL_iue,hgSL_jle,hgSL_jue,hgSL_kle,hgSL_kue,&
                         hgSR_ile,hgSR_iue,hgSR_jle,hgSR_jue,hgSR_kle,hgSR_kue)
        call Timers_stop("SolveBlock")        
        do k = 1, NZB
           do j = 1, NYB
              do i = 1, NXB
                 unk(gr_iSoln,i+NGUARD,j+K2D*NGUARD,k+K3D*NGUARD,b) =  soln(i,j,k)
                 work(i+NGUARD,j+K2D*NGUARD,k+K3D*NGUARD,b,1)       =  soln(i,j,k)
              enddo
           enddo
        enddo
        
     endif
  enddo
  if (nblockssolved == 0) then
     call Timers_start("fft")   !trick to keep timers structure on different procs the same - KW
     call Timers_stop("fft")
  end if

  call Timers_stop(level_timer_label)
  return
end subroutine gr_hgSolveLevel
