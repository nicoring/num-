module matrix
    implicit none
    public dot_product

  contains

    subroutine dot_product(A, B, C)
      real, intent(in) :: A(:, :)
      real, intent(in) :: B(:, :)
      real, intent(out), allocatable :: C(:, :)

      integer :: row
      integer :: col
      integer :: inner

      allocate(C(size(A, 1), size(B, 2)))

      do concurrent (row=1:size(A,1), col=1:size(B, 2), inner=1:size(A,2))
        C(row, col) = C(row, col) + A(row, inner) * B(inner, col)
      end do
  
    end subroutine dot_product
  
  end module matrix


program test_matrix
    use matrix
    implicit none

    real, allocatable :: A(:, :)
    real, allocatable :: B(:, :)
    real, allocatable :: C(:, :)
    
    allocate(A(100000, 1000))
    allocate(B(1000, 1000))

    A(:, :) = 0.0
    B(:, :) = 0.0

    call dot_product(A, B, C)
end program test_matrix