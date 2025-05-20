FC = gfortran

#FFLAGS = -g -fdefault-real-8 -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace \
#-ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan

FFLAGS = -g -fdefault-real-8 -fimplicit-none

LDFLAGS = 

SRCS = master.f90 mod_vgrid.f90 mod_consts.f90 mod_hgrid.f90 thermo_functions.f90 compute_t_qv_qc.f90 \
eddy_coefficients.f90 time_evolution.f90 pressure_gradx.f90 pressure_grady.f90 \
surface_fluxes.f90 advech_hq.f90 advech_uv.f90 diff_horiz_hq.f90 diff_horiz_uv.f90 initial_conditions_ideal.f90

OBJS = ${SRCS:.f90=.o}

all: main

main: $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJS) *.mod



