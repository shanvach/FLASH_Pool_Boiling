      subroutine init_1d (n1d,nvars,xzn,model_1d)
! 
! init_1d --
! 
!     read in the initial model from the input file and store it in the
!     model_1d array.  The initial model is stored in 1d arrays that are 
!     then used by init_block to fill the paramesh data structures
!
!    model_file -- the name of the file containing the model 1d atmosphere
!
!     for now, this assumes that the input file contains the x coordinate
!     and the 23 variables together
!  
      use dBase, ONLY :  dBasePropertyInteger

      use runtime_parameters

      implicit none

      integer i, j, n
      
      real dx

      integer n1d
      integer nvars
      
      real xzn(n1d), model_1d(n1d,nvars)

      integer op, kat

      real err

      integer fail, ierr
      parameter (fail = -1)

      integer          imax
      parameter        (imax=8192)


      integer itot

! storage for the data read in from the model file
      real  zdist(imax),ztemp(imax,nvars)

      character(len=80) :: model_file
      integer :: myPE, masterPE
      real :: xmin, xmax
!----------------------------------------------------------end of declarations
      call get_parm_from_context ("model_file", model_file)
      call get_parm_from_context ("xmin", xmin)
      call get_parm_from_context ("xmax", xmax)

      MyPE = dBasePropertyInteger("MyProcessor")
      MasterPE = dBasePropertyInteger("MasterProcessor")


! print out information
      if (MyPE .EQ. MasterPE) then 
         print *, 'reading in initial model'
         print *, 'nvars =',nvars
      endif
      

! Set up zones 
      dx = (xmax-xmin)/real(n1d)
      xzn(1) = xmin + dx/2.

      do i = 2,n1d
         xzn(i)  = xzn(i-1) + dx
      end do
      
! read in the model
      open(unit=2,file=trim(model_file),status='old')

      print *, 'openned input file: ', model_file

      itot = 0

! note- to really fix the following, the velocities
! should be set to zone centers
!
      do i=1,imax
         read(2,*,end=11) zdist(i), (ztemp(i,j),j=1,nvars)
         itot = itot + 1
      enddo

      print *, 'done reading initial model'
11    close(unit=2)

      print *, 'file = ', trim(model_file), '  number points =', itot

! define the order of the interpolation
      op  = 2
      kat = 1

! do the interpolations and fill the model_1d array with the model file
! information at the n1d points along the flame direction

      do i = 1, n1d

! find the location of this zone's center coordinate
         call hunt(zdist,itot,xzn(i),kat)
         kat = max(1, min(kat - op/2 + 1, itot - op + 1))

! loop over the nvars variables and interpolate from the input file into the
! model_1d array, using the xzn location.  Note that ztemp is stored such
! that the distance is the first variable, and thus adjacent in memory.
         do n = 1, nvars
            call polint(zdist(kat),ztemp(kat,n),op,xzn(i), & 
     &           model_1d(i,n),err)
         enddo

      enddo

! print out information
      if (MyPE .EQ. MasterPE) then 
         open(unit=98,file='test.out',status='unknown')
         open(unit=99,file='test2.out',status='unknown')
         do i = 1,n1d
         write(98,*)i
         write(99,101)xzn(i),model_1d(i,1:nvars)
         enddo
         close(98)
         close(99)
      endif
101   format(18(1E14.8,2x))

      return
      end




