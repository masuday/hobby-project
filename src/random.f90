module random

use,intrinsic :: iso_fortran_env,only:int32,int64,real32,real64
use random_number_engines, only: splitmix64
implicit none

interface
   subroutine random_number_engine(state,harvest)
      use,intrinsic :: iso_fortran_env, only:int64
      integer(int64),intent(inout),allocatable :: state(:)
      integer(int64),intent(out) :: harvest
   end subroutine random_number_engine

   subroutine random_number_state_initializer(state,seed)
      use,intrinsic :: iso_fortran_env, only:int64
      integer(int64),intent(inout) :: state(:)
      integer(int64),intent(in),optional :: seed(:)
   end subroutine random_number_state_initializer
end interface

! main object
type random_number_generator
   private
   integer :: state_size = 0
   integer(int64),allocatable :: state(:)
   procedure(random_number_engine),nopass,pointer :: &
      engine => random_number_engine_undef
   procedure(random_number_state_initializer),nopass,pointer :: &
      initializer => random_number_state_initializer_undef
end type random_number_generator

private
public :: random_number_generator, random_number_seed, random_real_uniform

contains

!
! Seeding
!
subroutine random_number_seed(rng, init_engine, size, put, put_scalar, get, auto)
   type(random_number_generator),intent(inout) :: rng   ! state structure
   character(len=*),intent(in),optional :: init_engine  ! engine name
   integer,intent(out),optional :: size                 ! same as the reference
   integer(int64),intent(in),optional :: put(:)         ! same as the reference
   integer(int64),intent(in),optional :: put_scalar     ! same as the reference
   logical,intent(in),optional :: auto                  ! automatic seeding
   integer(int64),intent(out),optional :: get(:)        ! same as the reference

   character(len=32) :: engine_name

   if(present(init_engine)) then
      select case(engine_name)
      case("default","xoshiro128**","xoshiro128plusplus")
         call random_generator_xoshiro128starstar(rng)
      case("mt19937_64","mt19937-64")
         call random_generator_mt19937_64(rng)
      case default
         stop "Unknown engine name."
      end select
   end if

   if(rng%state_size < 1) then
      stop "PRNG is not initialized. Call this subroutine with 'init_engine='."
   end if

   if(present(size)) then
      size = rng%state_size
   end if

   if(present(put)) then
      if(ubound(put,1)/=rng%state_size) then
         stop "wrong state size of 'put'"
      end if
      rng%state(1:rng%state_size) = put(1:rng%state_size)
      call rng%initializer(rng%state(1:rng%state_size),put(1:rng%state_size))
   end if

   if(present(put_scalar)) then
      rng%state(1:rng%state_size) = put_scalar
      call rng%initializer(rng%state(1:rng%state_size),[put_scalar])
   end if

   if(present(auto)) then
      if(auto) then
         call put_serial_date_time(rng%state(1:rng%state_size))
         call rng%initializer(rng%state(1:rng%state_size))
      end if
   end if

   if(present(get)) then
      if(ubound(get,1)/=rng%state_size) then
         stop "wrong state size of `get`"
      end if
      get(1:rng%state_size) = rng%state(1:rng%state_size)
   end if
end subroutine random_number_seed

!
! Assignment of an engine to the object
!
subroutine random_generator_xoshiro128starstar(rng)
   type(random_number_generator),intent(inout) :: rng
   rng%state_size = 2
   rng%state = [1,2]
   rng%engine => random_number_engine_xoshiro128starstar
   rng%initializer => random_number_state_initializer_splitmix64
   call rng%initializer(rng%state(1:rng%state_size))
end subroutine random_generator_xoshiro128starstar

subroutine random_generator_mt19937_64(rng)
   type(random_number_generator),intent(inout) :: rng
   ! keeping mti as the 313th element in the state
   rng%state_size = 312
   allocate(rng%state(1:313))
   rng%state(1:312) = 123456789_int64
   rng%state(1:313) = 313
   rng%engine => random_number_engine_xoshiro128starstar
   rng%initializer => random_number_state_initializer_mt19937_64
   call rng%initializer(rng%state(1:rng%state_size))
end subroutine random_generator_mt19937_64

subroutine random_number_engine_undef(state,harvest)
   integer(int64),intent(inout),allocatable :: state(:)
   integer(int64),intent(out) :: harvest
   stop "undefined engine"
end subroutine random_number_engine_undef

!
! Initializer
!
subroutine random_number_state_initializer_undef(state,seed)
   use,intrinsic :: iso_fortran_env, only:int64
   integer(int64),intent(inout) :: state(:)
   integer(int64),intent(in),optional :: seed(:)
   return
end subroutine random_number_state_initializer_undef

subroutine random_number_state_initializer_splitmix64(state,seed)
   use,intrinsic :: iso_fortran_env, only:int64
   integer(int64),intent(inout) :: state(:)
   integer(int64),intent(in),optional :: seed(:)
   integer :: i,n
   integer(real64) :: x,z,seed_size
   n = size(state)
   do i=1,n
      x = state(i)
      call splitmix64(x,z)
      state(i) = z
   end do
   return
end subroutine random_number_state_initializer_splitmix64

subroutine random_number_state_initializer_mt19937_64(state,seed)
   use,intrinsic :: iso_fortran_env, only:int64
   integer(int64),intent(inout) :: state(:)
   integer(int64),intent(in),optional :: seed(:)
   integer(int64) :: i,mti
   mti = 313
   call mt_init_by_array64(state(1:312),mti,seed)
   state(313) = mti
end subroutine random_number_state_initializer_mt19937_64

!
! Wrapper of random bit generator
!
subroutine random_number_engine_xoshiro128starstar(seed,rand)
   integer(int64),intent(inout),allocatable :: seed(:)
   integer(int64),intent(out) :: rand
   integer(int64) :: s0,s1
   s0 = seed(1)
   s1 = seed(2)
   rand = s0*5_int64
   rand = ior(ishft(rand,7),ishft(rand,-(64-7))) * 9_int64
   s1 = ieor(s1,s0)
   s0 = ieor(ieor(ior(ishft(s0,24),ishft(s0,-(64-24))),s1),ishft(s1,-16))
   s1 = ior(ishft(s1,37),ishft(s1,-(64-37)))
   seed(1) = s0
   seed(2) = s1
end subroutine random_number_engine_xoshiro128starstar

subroutine random_number_engine_mt19937_64(state,harvest)
   integer(int64),intent(inout),allocatable :: state(:)
   integer(int64),intent(out) :: harvest
   integer(int64) :: mti
   mti = state(313)
   call mt_genrand64_int64(state,mti,harvest)
   state(313) = mti
end subroutine random_number_engine_mt19937_64

!
! Uniform real PRNGs in [0,1)
!
subroutine random_real_uniform(rng,harvest)
   type(random_number_generator),intent(inout) :: rng
   real(real64),intent(out) :: harvest
   real(real64),parameter :: max_int64_double = 1/(real(huge(1_int64),kind=real64)-1)
   integer(int64) :: int_rand
   call rng%engine(rng%state,int_rand)
   harvest = abs(int_rand) * max_int64_double
end subroutine random_real_uniform

!
! Utilities
!
! taken from https://github.com/certik/hfsolver/blob/b4c50c1979fb7e468b1852b144ba756f5a51788d/src/utils.f90#L402
subroutine put_serial_date_time(seed)
   integer(int64),intent(inout) :: seed(:)
   integer :: i, n, clock
   n = size(seed)
   call system_clock(count=clock)
   seed = seed + clock + 37 * [(i - 1, i = 1, n)]
end subroutine put_serial_date_time

end module random


program main
   use,intrinsic :: iso_fortran_env,only:int32,int64,real64
   use random
   implicit none
   integer :: i
   real(real64) :: x
   type(random_number_generator) :: rng

   call random_number_seed(rng, "default", auto=.true.)
   do i=1,20
      call random_real_uniform(rng,x)
      print *,x
   end do
end program main
