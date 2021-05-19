############################################################
#                         MAKEFILE                         #
############################################################

# GNU Compiler
CMPF  = gfortran -c -cpp -I. -I/usr/local/include -g #-ffpe-trap=zero,invalid,overflow      # -fopenmp
LNK   = gfortran -llapack -llis-2.0-30 -L/students/st_19_20/vallisa/local/lib #-L/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/lapack/3.6.0/lib64 -lblas -L/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/lapack/3.6.0/lib64 #-libgfortran#-lc --entry main#gfortran  -L/usr/lib/x86_64-linux-gnu/lapack -llapack -L/usr/lib/x86_64-linux-gnu -llapacke  #-lgfortran# -L/usr/lib/x86_64-linux-gnu/lapack -llapack -L/usr/lib/x86_64-linux-gnu/blas -lblas  #-llapacke# -fopenmp
OPTF = -O0 -Wall -fimplicit-none -fbacktrace #-pg #-fopenmp

# Objects: list of all objects *.o
# OBJS = mpi_common.o  global.o  screen.o  tools.o  initialization.o  timecycle.o  EM_fields.o  grid_and_partition.o  particle.o  collisions.o  postprocess.o  mt19937.o

OBJS = time_integrators.o  numerical_fluxes.o  pde.o  pde_euler.o  pde_mom5.o  pde_mom14.o  pde_mom14_axi.o  various.o  EM_fields.o  grid.o  initialization.o  poisson_solver.o pde_euler_multi_fluid.o implicit_electrons_utilities.o rates_box.o global.o pde_euler_electrons.f90 pde_euler_ions.o

# Executable generation by the linker
formica.exe: formica.o $(OBJS)
	$(LNK) $(OPTF) formica.o $(OBJS) \
	            -o formica.exe

# Objects generation
formica.o: formica.f90 $(OBJS)
	$(CMPF) $(OPTF) formica.f90

time_integrators.o:  time_integrators.f90  EM_fields.o  numerical_fluxes.o  global.o poisson_solver.o implicit_electrons_utilities.o
	$(CMPF) $(OPTF) time_integrators.f90

numerical_fluxes.o:  numerical_fluxes.f90  pde.o  grid.o implicit_electrons_utilities.o global.o
	$(CMPF) $(OPTF) numerical_fluxes.f90

pde.o:  pde.f90 global.o pde_euler.o pde_mom5.o pde_mom14.o pde_mom14_axi.o pde_euler_multi_fluid.o
	$(CMPF) $(OPTF) pde.f90

pde_mom5.o:  pde_mom5.f90 global.o  various.o
	$(CMPF) $(OPTF) pde_mom5.f90

pde_mom14.o:  pde_mom14.f90 global.o  various.o
	$(CMPF) $(OPTF) pde_mom14.f90

pde_mom14_axi.o:  pde_mom14_axi.f90 global.o  various.o
	$(CMPF) $(OPTF) pde_mom14_axi.f90

pde_euler.o:  pde_euler.f90 global.o pde_euler_ions.o pde_euler_electrons.o
	$(CMPF) $(OPTF) pde_euler.f90

pde_euler_multi_fluid.o: pde_euler_multi_fluid.f90 global.o rates_box.o
	$(CMPF) $(OPTF) pde_euler_multi_fluid.f90

pde_euler_ions.o: pde_euler_ions.f90 global.o rates_box.o
	$(CMPF) $(OPTF) pde_euler_ions.f90

pde_euler_electrons.o: pde_euler_electrons.f90 global.o rates_box.o
	$(CMPF) $(OPTF) pde_euler_electrons.f90

EM_fields.o: EM_fields.f90 global.o
			$(CMPF) $(OPTF) EM_fields.f90

various.o: various.f90 EM_fields.o global.o
	$(CMPF) $(OPTF) various.f90

grid.o: grid.f90  global.o
	$(CMPF) $(OTPF) grid.f90

initialization.o: initialization.f90 global.o implicit_electrons_utilities.o rates_box.o EM_fields.o
	$(CMPF) $(OPTF) initialization.f90

poisson_solver.o: poisson_solver.f90 global.o EM_fields.o
	$(CMPF) $(OPTF) poisson_solver.f90

rates_box.o: rates_box.f90 global.o
	$(CMPF) $(OPTF) rates_box.f90

implicit_electrons_utilities.o: implicit_electrons_utilities.f90 pde_euler.o global.o pde_euler_electrons.o
	$(CMPF) $(OPTF) implicit_electrons_utilities.f90

global.o: global.f90
	$(CMPF) $(OTPF) global.f90

# Cleaning command
clean:
	@echo cleaning objects, modules and executables
	rm  -f  *.o  *.mod  *.exe  *~

cleanoutput:
	@echo cleaning output and dump files
	rm  -f  dumps/*
