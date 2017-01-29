!*******************************************************************************

! Routine:      readSdTable

! Description:  Reads in one of the Sutherland & Dopita (1993) files containing
!               plasma cooling curve data.  These files are available at
!               http://msowww.anu.edu.au/~ralph/data/cool/.

! Reference:    Sutherland, R. S., and Dopita, M. A. 1993, ApJS, 88, 253

! Parameters:   file              Name of file to read.
!               logT(:)           1D array returning log temperature in K.
!               ne(:)             Array returning electron number density
!                                   in cm^-3
!               nH(:)             Array returning hydrogen number density
!                                   in cm^-3
!               nt(:)             Array returning total ion number
!                                   density in cm^-3
!               logLnet(:)        Array returning cooling function Lambda
!                                   in log(erg cm^-3 s^-1)
!               logLnorm(:)       Array returning log (Lambda/(ne nt))
!               logU(:)           Array returning log(3/2(ne+nt)*kT)
!               logtau(:)         Array returning norm cooling timescale
!                                   log( U/(Lambda*(ne+nt)) )
!               P12               Array returning pressure in erg cm^-3
!                                   * 1.E12
!               rho24             Array returning density in g cm^-3
!                                   * 1.E24
!               Ci                Array returning isothermal sound speed
!                                   in km s^-1
!               mubar             Array returning mean molecular weight
!                                   in g * 1.E24
!               Nmax              Size of 1D arrays.
!               N                 Integer returning number of data rows found
!                                   in the file.  If > 0, this will be the
!                                   size of the arrays allocated for each of
!                                   the above pointers.  If < 0, there was a
!                                   problem and no data were read.


subroutine readSdTable (file, logT, ne, nH, nt, logLnet, logLnorm, &
                          logU, logtau, P12, rho24, Ci, mubar, Nmax, N)

!===============================================================================

  implicit none
  
  character(len=*), intent(IN) :: file
  integer, intent(IN)   :: Nmax
  integer, intent(OUT)  :: N

  real, intent(INOUT), dimension(Nmax) :: logT, ne, nH, nt, logLnet, logLnorm, &
                                          logU, logtau, P12, rho24, Ci, mubar
  
  real    :: row(12)
  integer :: i, status
  
  !===============================================================================
  
  ! First, determine whether the file can be opened and has the right format.
  ! Count the data rows in the file.
  
  open (1, file=file, iostat=status)
  if (status /= 0) then
     N = -1
     return
  endif
  
  do i = 1, 4    ! first 4 lines contain header info
     read (1,*)
  enddo
  
  N = 0
  
  do while (.true.)
     read (1,*,iostat=status) row
     if (status == -1) exit
     if (status > 0) then
        N = -1
        close (1)
        return
     endif
     N = N + 1
  enddo
  
  close (1)
  
  if ((N == 0) .or. (N > Nmax)) then
     N = -1
     return
  endif
  
  ! Next, free up any space that might have been allocated previously.
  
  !if (associated(logT))     deallocate (logT)
  !if (associated(ne))       deallocate (ne)
  !if (associated(nH))       deallocate (nH)
  !if (associated(nt))       deallocate (nt)
  !if (associated(logLnet))  deallocate (logLnet)
  !if (associated(logLnorm)) deallocate (logLnorm)
  !if (associated(logU))     deallocate (logU)
  !if (associated(logtau))   deallocate (logtau)
  !if (associated(P12))      deallocate (P12)
  !if (associated(rho24))    deallocate (rho24)
  !if (associated(Ci))       deallocate (Ci)
  !if (associated(mubar))    deallocate (mubar)
  
  ! Now allocate space for the arrays.
  
  !allocate (logT(N))
  !allocate (ne(N))
  !allocate (nH(N))
  !allocate (nt(N))
  !allocate (logLnet(N))
  !allocate (logLnorm(N))
  !allocate (logU(N))
  !allocate (logtau(N))
  !allocate (P12(N))
  !allocate (rho24(N))
  !allocate (Ci(N))
  !allocate (mubar(N))
  
  !===============================================================================
  
  ! Open the file again and read the data.
  
  open (1, file=file, iostat=status)
  
  do i = 1, 4    ! first 4 lines contain header info
     read (1,*)
  enddo
  
  do i = 1, N
     read (1,*) logT(i), ne(i), nH(i), nt(i), logLnet(i), logLnorm(i), &
          logU(i), logtau(i), P12(i), rho24(i), Ci(i), mubar(i)
  enddo
  
  close (1)
  
  !===============================================================================
  
  return
end subroutine readSdTable
