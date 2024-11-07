# GLAS Buildsystem Database
# This file is generated automatically by "install.py" at 2023-03-04

blas = glas_buildsystem.library()
blas.library = '-L/usr/lib -lblas'

lapack = glas_buildsystem.library()
lapack.library = '-L/usr/lib -llapack'
lapack += blas

arpack = glas_buildsystem.library()
arpack.library = '-L/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//ARPACK/ -larpack -L/usr/local/gfortran/lib -lgfortran'
arpack += lapack

boost = glas_buildsystem.library()
boost.include_path = ["/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//boost_1_59_0",]
boost.library = '-L/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//boost_1_59_0/libs -lboost_1_59_0'

boost_numeric_bindings = glas_buildsystem.library()
boost_numeric_bindings.include_path = ["/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software/",]
boost_numeric_bindings += boost

arpack_bindings = glas_buildsystem.library()
arpack_bindings.include_path = ["/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//arpack_bindings",]
arpack_bindings.library = '-L/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//arpack_bindings/lib -larpack_bindings'
arpack_bindings += arpack
arpack_bindings += boost_numeric_bindings

mumps = glas_buildsystem.library()
mumps.include_path = ["/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//MUMPS_5.0.1/include","/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//MUMPS_5.0.1/libseq","/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//MUMPS_5.0.1/pord/include","/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//MUMPS_5.0.1/src",]
mumps.library = '-L/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//MUMPS_5.0.1/libseq -L/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//MUMPS_5.0.1/lib -lcmumps -ldmumps -lmumps_common -lzmumps -lpord -lmpiseq -L/usr/local/gfortran/lib -lgfortran'
mumps += lapack

glas2 = glas_buildsystem.library()
glas2.include_path = ["/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software/glas2",]
glas2.library = '-L/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software/glas2/libs/ -lglas2'
glas2.defines = '-DGLAS_COMPLEX'
glas2 += boost

itsol = glas_buildsystem.library()
itsol.include_path = ["/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//ITSOL_2/INC",]
itsol.library = '-L/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//ITSOL_2 -litsol'

itsol_pp = glas_buildsystem.library()
itsol_pp.include_path = ["/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//itsol++",]
itsol_pp.library = '-L/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//itsol++/lib/ -litsol_pp'
itsol_pp += itsol
itsol_pp += boost_numeric_bindings

