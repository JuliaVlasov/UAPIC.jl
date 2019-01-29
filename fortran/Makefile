SRCS =	compute_rho.f90 interpolation.f90 meshfields.f90 particles.f90 \
	    poisson.f90 ua_steps.f90 ua_type.f90 

OBJS =	compute_rho.o interpolation.o meshfields.o particles.o poisson.o \
	    ua_steps.o ua_type.o

LIBS =	-lfftw3

F90 = gfortran
F90FLAGS = -I/usr/local/include -fstack-arrays -O0 -g -fbacktrace -pedantic -Wuninitialized -fcheck=array-temps,bounds,do,pointer,recursion -ffpe-trap=invalid,zero,overflow
F90FLAGS = -I/usr/local/include -O3
LDFLAGS = $(F90FLAGS)  

all: bupdate efd

bupdate: bupdate.f90 $(OBJS)
	$(F90) $(LDFLAGS) $< -o $@ $(OBJS) $(LIBS)

efd: efd.f90 $(OBJS)
	$(F90) $(LDFLAGS) $< -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.mod.o:
	$(F90) $(F90FLAGS) -c $*.f90

meshfields.o: meshfields.f90
particles.o: particles.f90 meshfields.o
compute_rho.o: compute_rho.f90 meshfields.o particles.o ua_type.o
interpolation.o: interpolation.f90 meshfields.o particles.o compute_rho.o
poisson.o: poisson.f90 meshfields.o 
ua_steps.o: ua_steps.f90 compute_rho.o meshfields.o particles.o ua_type.o
ua_type.o: ua_type.f90 
