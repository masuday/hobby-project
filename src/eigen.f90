!
! Eigenvector of a real, symmetric matrix.
! Method by Denton et al. (https://arxiv.org/abs/1908.03795)
! Program by Yutaka Masuda
!
! This program is in public domain.
!
program eigen
   use iso_fortran_env

   implicit none
   character(len=16) :: argv
   integer :: n,i,j,k,lwork,info
   real(real64) :: lhs,rhs,ev
   real(real64),allocatable :: A(:,:),M(:,:),evalA(:),evalM(:),evecA(:),work(:)
   external :: DSYEV,ILAENV

   ! n = size; j = excluded row/column
   if(command_argument_count()<2) then
      stop 'usage: eigen n j'
   end if
   call get_command_argument(1,argv)
   read(argv,*,iostat=info) n
   if(info/=0 .or. n<1) stop 'invalid n = size(A)'
   call get_command_argument(2,argv)
   read(argv,*,iostat=info) j
   if(info/=0 .or. j<1 .or. j>n) stop 'invalid j = excluded row/column <= n'

   ! preparation for matrix
   allocate(A(n,n),M(n-1,n-1),evalA(n),evalM(n-1),evecA(n))
   call random_number(A)
   A = A + transpose(A)
   call make_submatrix(A,M,j)

   ! for LAPACK
   lwork = get_optimal_worksize(n)
   allocate(work(lwork))

   ! eigenvalues/vectors of full matrix (A)
   call DSYEV('V','L',n,A,n,evalA,work,lwork,info)
   if(info/=0) stop 'failed in DSYEV for A'

   ! eigenvalues of submatrix (M)
   call DSYEV('N','L',n-1,M,n-1,evalM,work,lwork,info)
   if(info/=0) stop 'failed in DSYEV for M'

   ! new algorithm
   call get_eigenvector(evalA,evalM,evecA)
   if(n<=10) then
      print '(A,*(ES16.8))','traditional',A(j,:)**2
      print '(A,*(ES16.8))','new method ',evecA
   end if
   print '(A,*(ES16.8))','max|diff|  ',maxval(abs(A(j,:)**2 - evecA))

contains

subroutine make_submatrix(A,M,j)
   real(real64),intent(in) :: A(:,:)
   real(real64),intent(out) :: M(:,:)
   integer,intent(in) :: j
   integer,allocatable :: idx(:)
   integer :: i,n
   n = size(A,1)
   if(size(A,1)/=size(A,2)) stop 'A not square'
   if(size(M,1)/=size(M,2)) stop 'M not square'
   if(size(M,1)/=n-1) stop 'size(M,1) /= size(A,1)-1'
   if(j<1 .or. j>n) stop 'k out of range'
   idx = [integer::]
   do i=1,n
      if(i/=j) idx=[idx,i]
   end do
   M(1:n-1,1:n-1) = A(idx,idx)
   deallocate(idx)
end subroutine make_submatrix

subroutine get_eigenvector_naive(evalA,evalM,evecA)
   real(real64),intent(in) :: evalA(:),evalM(:)
   real(real64),intent(inout) :: evecA(:)
   integer :: i,k,n
   real(real64) :: lhs,rhs,ev
   n = size(evalA)
   do i=1,n
      lhs = 1.0d0
      do k=1,n
         if(k/=i) then
            lhs = lhs * (evalA(i)-evalA(k))
         end if
      end do
      rhs = 1.0d0
      do k=1,n-1
         rhs = rhs * (evalA(i)-evalM(k))
      end do
      evecA(i) = rhs/lhs
   end do
end subroutine get_eigenvector_naive

subroutine get_eigenvector(evalA,evalM,evecA)
   real(real64),intent(in) :: evalA(:),evalM(:)
   real(real64),intent(inout) :: evecA(:)
   integer :: i,k,n
   real(real64) :: lhs,rhs,ev,d,s,t
   n = size(evalA)
   do i=1,n
      lhs = 0.0d0
      s = 1.0d0
      do k=1,n
         if(k/=i) then
            d = evalA(i) - evalA(k)
            s = s * sign(1.0d0,d)
            lhs = lhs + log(abs(d))
         end if
      end do
      rhs = 0.0d0
      t = 1.0d0
      do k=1,n-1
         d = evalA(i) - evalM(k)
         t = t * sign(1.0d0,d)
         rhs = rhs + log(abs(d))
      end do
      evecA(i) = s*t*exp(rhs-lhs)
   end do
end subroutine get_eigenvector

function get_optimal_worksize(n) result(lwork)
   integer,intent(in) :: n
   integer :: nb,n1,n2,n3,n4,lwork
   integer :: ILAENV
   external :: ILAENV
   nb = ILAENV(1,'DSYTRD','',n1,n2,n3,n4)
   lwork = (nb+2)*n
end function get_optimal_worksize

subroutine print_matrix(A,msg)
   real(real64),intent(in) :: A(:,:)
   character(len=*),intent(in),optional :: msg
   integer :: i,j
   if(present(msg)) print '(A)',trim(msg)
   do i=1,size(A,1)
      print '(*(ES16.8))',(A(i,j),j=1,size(A,2))
   end do
end subroutine print_matrix

end program eigen
