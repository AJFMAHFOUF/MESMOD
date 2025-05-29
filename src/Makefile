FC = gfortran

#FFLAGS = -g -fdefault-real-8 -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace \
#-ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan

FFLAGS = -g -fdefault-real-8 -fimplicit-none

LDFLAGS = 

SRCS = master.f90 mod_vgrid.f90 mod_consts.f90 mod_hgrid.f90 advecv_h.f90 advecv_u.f90 \
compute_geopotential.f90 time_evolution.f90 pressure_gradx.f90  \
advech_h.f90 advech_u.f90 diff_horiz_h.f90 diff_horiz_u.f90 

OBJS = ${SRCS:.f90=.o}

all: main

main: $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJS) *.mod



