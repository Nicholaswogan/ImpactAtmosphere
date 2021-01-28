FCFLAGS="-O3"

gfortran -c src/minpack/dpmpar.f $FCFLAGS
gfortran -c src/minpack/enorm.f $FCFLAGS
gfortran -c src/minpack/fdjac2.f $FCFLAGS
gfortran -c src/minpack/lmdif.f $FCFLAGS
gfortran -c src/minpack/lmdif1.f $FCFLAGS
gfortran -c src/minpack/lmdif2.f $FCFLAGS
gfortran -c src/minpack/lmpar.f $FCFLAGS
gfortran -c src/minpack/qrfac.f $FCFLAGS
gfortran -c src/minpack/qrsolv.f $FCFLAGS

gfortran -c src/diffusion.f90 $FCFLAGS
gfortran -c src/EvolveAtmFort.f90 $FCFLAGS

# gfortran EvolveAtmFort.o diffusion.o dpmpar.o enorm.o fdjac2.o lmdif.o lmdif1.o lmdif2.o lmpar.o qrfac.o qrsolv.o -o test.run $FCFLAGS -llapack


python -m numpy.f2py -c src/diffusion.f90 src/EvolveAtmFort.f90 src/minpack/dpmpar.f src/minpack/enorm.f \
                        src/minpack/fdjac2.f src/minpack/lmdif.f src/minpack/lmdif1.f \
                        src/minpack/lmdif2.f src/minpack/lmpar.f src/minpack/qrfac.f src/minpack/qrsolv.f \
-m EvolveAtmFort \
--opt="-O3" \
-llapack \
only: setup rhs rhs_verbose diffuse


rm *.mod
rm *.o
