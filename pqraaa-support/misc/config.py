# C++ compiler
cxx = glas_buildsystem.new_compiler( 'C++', ['cpp'] )
cxx.debug = '-g'
cxx.release = '-DNDEBUG -O3 -funroll-loops'
cxx.compile = 'c++ -c -std=c++17 -Wall '#-DGLAS_OPENMP ' #-DBIND_FORTRAN_LOWERCASE '
cxx.link = 'c++ -o'
cxx.library = '-lstdc++ -lpthread'
cxx.ar = 'ar -r'

# Fortran compiler (only needed when LAPACK and BLAS are not installed)
fortran = glas_buildsystem.new_compiler('Fortran', ['f','F','f90'])
fortran.debug = '-g'
fortran.release = '-O3 -funroll-loops'
fortran.compile = 'gfortran -c -Wall '
fortran.library = '-L/usr/local/gfortran/lib -lgfortran'
# Fortran names get an underscore for linking to C
fortran.c_name_ext = 'Add_'

c = glas_buildsystem.new_compiler('C', ['c'])
c.debug = '-g'
c.release = '-O3 -funroll-loops'
c.compile = 'gcc -c -Wall '

# BLAS (Remove or comment if not installed)
blas = glas_buildsystem.library()
blas.library = '-L/usr/lib -lblas'

# LAPACK (depends on BLAS) (Remove or comment if not installed)
lapack = glas_buildsystem.library()
lapack.library = '-L/usr/lib -llapack'

# Add openmp if multithreading is going to be used.
# Openmp is needed to use the openmp backend in glas2
openmp = glas_buildsystem.library()
openmp.flags = '-pthread -fopenmp'
openmp.library = '-pthread -fopenmp'

# Boost (if available)
#boost = glas_buildsystem.library()
#boost.include = ['/usr/local/include']
#boost.library = '-L/usr/lib -lboost'
