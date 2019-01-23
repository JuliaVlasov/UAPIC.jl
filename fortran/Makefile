PROG = uapic_2d	

SRCS =	$(wildcard *.F90)

OBJS =	$(SRCS:.F90=.o) 

LIBS =	-lfftw3

F90 = gfortran
F90FLAGS = -I/usr/local/include -g -C
LDFLAGS = 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .F90

.F90.o:
	$(F90) $(F90FLAGS) -c $<

.mod.o:
	$(F90) $(F90FLAGS) -c $*.F90

ua_type.o: particles.o meshfields.o
compute_rho.o: ua_type.o 
