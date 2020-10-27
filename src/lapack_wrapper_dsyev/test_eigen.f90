!
! A test program for eigendecomposition.
! using a custom wrapper for LAPACK.
!
program test_eigen
   use lapack_wrapper
   implicit none
   integer,parameter :: dp=selected_real_kind(15,300)
   integer,parameter :: n=5
   integer :: info
   real(dp) :: a(n,n),x(n,n),eigval(n)

   ! make a random symmetric matrix
   call random_number(a)
   a = transpose(a) + a
   call show(a, "original matrix:")

   ! The eigen subroutine will destroy the input matrix.
   ! Make a copy if you want to keep the original matrix.
   ! See lapack_wrapper.f90 for details.

   ! Get the eigenvalues only.
   ! The variable x will be destroyed.
   print '(/,a)','Eigenvalue decomposition: case 1:'
   x = a
   call dsyev_f95sub(x, eigval, "N", "L", info)
   if(info/=0) then
      call error_stop("Error in eigendecomposition.")
   end if
   print '(a)','Eigenvalues:'
   print '(*(f12.8))',eigval(1:n)

   ! Get both eigenvalues and eigenvectors.
   ! The variable x will be updated with a set of eigenvectors.
   print '(/,a)','Eigenvalue decomposition: case 2:'
   x = a
   call dsyev_f95sub(x, eigval, "V", "L", info)
   if(info/=0) then
      call error_stop("Error in eigendecomposition.")
   end if
   print '(a)','Eigenvalues:'
   print '(*(f12.8))',eigval(1:n)
   call show(x, "Eigenvectors:")

   ! Test if the decomposition is valid
   call test_eigendecomposition(a, eigval,x)

contains

subroutine test_eigendecomposition(a, eigval,eigvec)
   real(dp),intent(in) :: a(:,:),eigval(:),eigvec(:,:)
   real(dp) :: avgdiff,diag(size(eigval),size(eigval))
   integer :: i,n

   diag = 0.0
   n = size(eigval)
   do i=1,n
      diag(i,i) = eigval(i)
   end do

   ! Test if avg|A - V*D*V'| is close to 0.
   avgdiff = sum(abs(a - matmul(matmul(eigvec,diag),transpose(eigvec)))) / (n*n)
   if(avgdiff > 1e-8) then
      print '(a,es12.4)','Failed.',avgdiff
   else
      print '(a,es12.4)','Okay.  ',avgdiff
   end if
end subroutine test_eigendecomposition

subroutine show(a,msg)
   character(len=*),intent(in),optional :: msg
   real(dp),intent(in) :: a(:,:)
   integer :: i
   if(present(msg)) then
      print '(a)',msg
   end if
   do i=1,size(a,1)
      print '(*(f12.8))',a(i,:)
   end do
end subroutine show

subroutine error_stop(msg)
   character(len=*),intent(in),optional :: msg
   if(present(msg)) then
      print '(a)',msg
   end if
   stop
end subroutine error_stop

end program test_eigen
