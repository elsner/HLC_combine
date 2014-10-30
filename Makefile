
FC      = mpif90
#FC     += -O0 -openmp -traceback -warn all -check all
FC     += -O3 -openmp

INCDIR  = /path/to/Healpix/include
LIBDIR  = /path/to/libhealpix/libcfitsio

FFLAGS  = -I$(INCDIR)
LDFLAGS = -L$(LIBDIR) -lhealpix -lcfitsio

exec    = HLC_combine
src     = combine_main.f90
obj     = combine_modules.o combine_main.o


default: $(src) $(exec)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(exec): $(obj)
	$(FC) -o $(exec) $(obj) $(LDFLAGS)

clean:
	rm -f *~ *.o *.mod

cleanall: clean
	rm -f $(exec)
