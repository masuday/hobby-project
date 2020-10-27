module random_number_engines

use,intrinsic :: iso_fortran_env, only:int64

implicit none

! for Mersenne Twister 64
integer(int64),parameter,private :: NN = 312_int64
integer(int64),parameter,private :: MM = 156_int64
integer(int64),parameter,private :: MATRIX_A = -5403634167711393303_int64
integer(int64),parameter,private :: UM = -2147483648_int64
integer(int64),parameter,private :: LM = 2147483647_int64

private :: rotl
public :: &
    xoroshiro128plusplus, xoroshiro128starstar, &
    xoshiro256plusplus, xoshiro256starstar, &
    splitmix64, &
    mt_init_genrand64, mt_init_by_array64, mt_genrand64_int64

contains

!
! source:
! http://prng.di.unimi.it/xoroshiro128plusplus.c
! Written in 2019 by David Blackman and Sebastiano Vigna
! public domain
!
subroutine xoroshiro128plusplus(s,harvest)
    integer(int64),intent(inout) :: s(0:1)
    integer(int64),intent(out) :: harvest
    integer(int64) :: s0,s1
    s0 = s(0)
    s1 = s(1)
    harvest = rotl(s0+s1,17) + s0
    s1 = ieor(s1,s0)
    s(0) = ieor(ieor(rotl(s0,49),s1),ishft(s1,21))
    s(1) = rotl(s1,28)
end subroutine xoroshiro128plusplus

!
! source:
! http://prng.di.unimi.it/xoroshiro128starstar.c
! Written in 2018 by David Blackman and Sebastiano Vigna
! public domain
!
subroutine xoroshiro128starstar(s,harvest)
    integer(int64),intent(inout) :: s(0:1)
    integer(int64),intent(out) :: harvest
    integer(int64) :: s0,s1
    s0 = s(0)
    s1 = s(1)
    harvest = rotl(s0*5_int64,7) * 9_int64
    s1 = ieor(s1,s0)
    s(0) = ieor(ieor(rotl(s0,24),s1),ishft(s1,16))
    s(1) = rotl(s1,37)
end subroutine xoroshiro128starstar

!
! source:
! http://prng.di.unimi.it/xoroshiro256plusplus.c
! Written in 2019 by David Blackman and Sebastiano Vigna
! public domain
!
subroutine xoshiro256plusplus(s,harvest)
    integer(int64),intent(inout) :: s(0:3)
    integer(int64),intent(out) :: harvest
    integer(int64) :: t
    harvest = rotl(s(0)+s(3),23) + s(0)
    t = ishft(s(1),17)
    s(2) = ieor(s(2),s(0))
    s(3) = ieor(s(3),s(1))
    s(1) = ieor(s(1),s(2))
    s(0) = ieor(s(0),s(3))
    s(2) = ieor(s(2),t)
    s(3) = rotl(s(3),45)
end subroutine xoshiro256plusplus

!
! source:
! http://prng.di.unimi.it/xoroshiro256starstar.c
! Written in 2018 by David Blackman and Sebastiano Vigna
! public domain
!
subroutine xoshiro256starstar(s,harvest)
    integer(int64),intent(inout) :: s(0:3)
    integer(int64),intent(out) :: harvest
    integer(int64) :: t
    harvest = rotl(s(1)*5_int64,7)
    t = ishft(s(1),17)
    s(2) = ieor(s(2),s(0))
    s(3) = ieor(s(3),s(1))
    s(1) = ieor(s(1),s(2))
    s(0) = ieor(s(0),s(3))
    s(2) = ieor(s(2),t)
    s(3) = rotl(s(3),45)
end subroutine xoshiro256starstar

function rotl(x,k) result(z)
    integer(int64),intent(in) :: x
    integer,intent(in) :: k
    integer(int64) :: z
    z = ior(ishft(x,k),ishft(x,-(64-k)))
end function rotl

!
! source:
! http://prng.di.unimi.it/splitmix64.c
! Written in 2015 by Sebastiano Vigna
! public domain
!
subroutine splitmix64(x,z)
    integer(int64),intent(inout) :: x
    integer(int64),intent(inout) :: z
    x = x - 7046029254386353131_int64
    z = x
    z = ieor(z,ishft(z,-30)) * (-4658895280553007687_int64)
    z = ieor(z,ishft(z,-27)) * (-7723592293110705685_int64)
    z = ieor(z,ishft(z,-31))
end subroutine splitmix64

!
! source:
! http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c
!   A C-program for MT19937-64 (2004/9/29 version).
!   Coded by Takuji Nishimura and Makoto Matsumoto.
! The 3-clause BSD license
!
subroutine mt_init_genrand64(mt,mti,seed)
    integer(int64),intent(inout) :: mt(0:)   ! state
    integer(int64),intent(inout) :: mti
    integer(int64),intent(in) :: seed
    integer(int64) :: i
    mt(0) = seed
    do i=1,NN-1
        mt(i) = (6364136223846793005_int64 * ieor(mt(i-1), ishft(mt(i-1), -62)) + i)
    end do
    mti = i
end subroutine mt_init_genrand64

subroutine mt_init_by_array64(mt,mti,init_key)
    integer(int64),intent(inout) :: mt(0:)   ! state
    integer(int64),intent(inout) :: mti
    integer(int64),intent(in) :: init_key(0:)
    integer(int64),parameter :: v1=3935559000370003845_int64
    integer(int64),parameter :: v2=2862933555777941757_int64
    integer(int64) :: size_key,i,j,k,kx

    i = 1
    j = 0
    size_key = size(init_key)
    call mt_init_genrand64(mt,mti,19650218_int64)
    k = max(NN,size_key)
    do kx=k,1,-1
       mt(i) = ieor(mt(i),v1*ieor(mt(i-1),ishft(mt(i-1),-62))) &
             + init_key(j) + j
       i = i + 1
       j = j + 1
       if(i>=NN) then
           mt(0) = mt(NN-1)
           i = 1
       end if
       if(j>=size_key) then
           j = 0
       end if
    end do
    do kx=nn-1,1,-1
        mt(i) = ieor(mt(i),v2*ieor(mt(i-1),ishft(mt(i-1),-62))) &
              - i
        i = i + 1
        if(i>=NN) then
            mt(0) = mt(NN-1)
            i = 1
        end if
    end do
    mt(0) = ishft(1_int64,63)
end subroutine mt_init_by_array64

subroutine mt_genrand64_int64(mt,mti,harvest)
    integer(int64),intent(inout) :: mt(0:)   ! state
    integer(int64),intent(inout) :: mti
    integer(int64),intent(out) :: harvest
    integer(int64),parameter :: mag01(0:1) = [0_int64, MATRIX_A]
    integer(int64) :: i,x
    if(mti>=NN) then
        if(mti==NN+1) then
            call mt_init_genrand64(mt,mti,5489_int64);
        end if
        do i=0,NN-MM-1
            x = ior(iand(mt(i),UM),iand(mt(i+1),LM))
            mt(i) = ieor(ieor(mt(i+mm),ishft(x,-1)), mag01(iand(x,1_int64)))
        end do
        do i=NN-MM,NN-1
            x = ior(iand(mt(i),UM),iand(mt(i+1),LM))
            mt(i) = ieor(ieor(mt(i+(MM-NN)),ishft(x,-1)), mag01(iand(x,1_int64)))
        end do
        x = ior(iand(mt(NN-1),UM),iand(mt(0),LM))
        mt(NN-1) = ieor(ieor(mt(MM-1),ishft(x,-1)), mag01(iand(x,1_int64)))
        mti = 0
    end if
    x = mt(mti)
    mti = mti + 1
    x = ieor(x,iand(ishft(x,-29),6148914691236517205_int64))
    x = ieor(x,iand(ishft(x, 17),8202884508482404352_int64))
    x = ieor(x,iand(ishft(x, 37),  -2270628950310912_int64))
    x = ieor(x,ishft(x, -43))
    harvest = x
end subroutine mt_genrand64_int64

end module random_number_engines


!program test
!   use,intrinsic :: iso_fortran_env, only:int64
!   use random_number_engines
!   implicit none
!   integer(int64) :: init_key(4)=[ 74565,  144470,  214375,  284280]
!   integer(int64) :: seed(1:313),harvest,i,j,mti
!   mti = 313
!   call mt_init_by_array64(seed,seed(313),init_key)
!   do i=1,1000
!      call mt_genrand64_int64(seed,seed(313),harvest)
!      write(*,"(*(z20))",advance="no") harvest
!      if(mod(i,5)==0) write(*,*)
!   end do
!end program test
