subroutine gr_pfftTriDiag (lower, main, upper, rhs, x, length, ierr)

  implicit none

  real, dimension(length), intent(in)  :: lower, main, upper, rhs
  real, dimension(length), intent(out) :: x
  integer, intent(in)     :: length
  integer, intent(out)    :: ierr
  real, dimension(length) :: work
  real :: const
  real, parameter :: eps = 1.0e-16
  integer :: i

  ! is the matrix singular
  ierr = 0

  ! forward elimination
  const = main(1)
  x(1) = rhs(1) / main(1)
  do i = 2, length
    work(i) = upper(i-1) / const
    const = main(i) - lower(i) * work(i)

    ! if singular prevent div-by-zero
    if (const == 0.0 ) then
      const = const + eps
      ierr = 1
    endif

    x(i) = (rhs(i) - lower(i) * x(i-1)) / const
  end do


  ! backward substitution
  do i = length-1, 1, -1
    x(i) = x(i) - work(i+1) * x(i+1)
  end do
    
  return
 
end subroutine gr_pfftTriDiag

subroutine gr_pfftCyclicTriDiag(lower, main, upper, alpha, beta, rhs, x, length, ierr)

  implicit none

  real, dimension(length), intent(in)  :: lower, main, upper, rhs
  real, dimension(length), intent(out) :: x
  real, intent(in)    :: alpha, beta
  integer, intent(in) :: length
  integer, intent(out)    :: ierr
  real, dimension(length) :: prime, u, z
  real    :: delta, fact
  integer :: i, err

  ! is the matrix singular
  ierr = 0

  !
  ! Solve none periodic system
  !  
  call gr_pfftTridiag(lower, main, upper, rhs, x, length, err)
  ierr = ierr + err
   
  !
  ! Sherman-Morrison formula applied 
  !        to the tridiagonal system 
  !
  
  ! setup main diagonal of modified system
  delta = -main(1)
  prime(1) = main(1) - delta
  prime(2:length-1) = main(2:length-1)
  prime(length) = main(length) - alpha * beta / delta

  ! solve the system -->  A' x = rhs
  call gr_pfftTridiag(lower, prime, upper, rhs, x, length, err)
  ierr = ierr + err

  ! form vector u
  u(1) = delta
  u(2:length-1) = 0.
  u(length) = alpha

  ! solve the system -->  A' z = u
  call gr_pfftTridiag(lower, prime, upper, u, z, length, err)
  ierr = ierr + err

  ! form the coefficient -->  v x / (1 + v z) 
  fact = (x(1) + beta * x(length) / delta) / (1.0 + z(1) + beta * z(length) / delta)

  ! get the solution to the cyclic system
  x(:) = x(:) - fact * z(:)
  
  return

end subroutine gr_pfftCyclicTriDiag
