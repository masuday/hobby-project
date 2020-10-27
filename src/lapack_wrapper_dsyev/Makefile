# For Intel Fortran
#FC = ifort
#FCFLAGS = -O2 -mkl

# For GFortran
FC = gfortran
FCFLAGS = -O2
LDFLAGS = -L. -I. -llapack -lrefblas

a.out: lapack_wrapper.o test_eigen.o
	$(FC) lapack_wrapper.o test_eigen.o -o $@ $(LDFLAGS)

lapack_wrapper.o: lapack_wrapper.f90
	$(FC) $(FCFLAGS) -c $< -o $@

test_eigen.o: test_eigen.f90 lapack_wrapper.o
	$(FC) $(FCFLAGS) -c $< -o $@
