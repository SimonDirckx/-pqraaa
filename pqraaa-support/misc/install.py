#! /usr/bin/python

import os
import os.path
import sys
import shutil

home = os.getcwd() + '/'
sys.path.insert(1, home + 'glas2/py')
library_file = home + 'libraries.py'

import glas_buildsystem_database
database = glas_buildsystem_database.database( library_file )

import glas_buildsystem
import glas_test
import dep_tools

# Build makefile include from config file
import glas_buildsystem
f = open('config.py','r')
exec( f.read() )
f.close()

if not 'cxx' in locals():
  print( "Help: I need the C++ compiler. Please add the following to config.py:")
  print ("cxx.debug = '-g' # or other options for debugging purposes")
  print ("cxx.release = '-DNDEBUG -O3 -funroll-loops -fvectorize' # or other options for optimization")
  print ("cxx.compile = 'c++ -c -std=c++11 -Wall ' # or other program for compiling")
  print ("cxx.link = 'c++ -o' # or other program for linking")
  print ("cxx.library = '-lstdc++' # or other libraries to link C++ programs with")

if not 'fortran' in locals():
  print ('Help: I need the Fortran (77+90) compiler. Please add the following to config.py:')
  print ("fortran = glas_buildsystem.new_compiler('Fortran', ['f','F'])")
  print ("fortran.debug = '-g' # or other options for debugging purposes")
  print ("fortran.release = '-O3 -funroll-loops' # or other options for optimization")
  print ("fortran.compile = 'gfortran -c -Wall ' # or other program for compiling")
  quit('Please add fortran compiler to config file')
  print ("cxx.ar = 'ar -r' # or other command for making library archives")

if not 'c' in locals():
  print ('Help: I need the C compiler. Please add the following to config.py:')
  print ("c = glas_buildsystem.new_compiler('C', ['c'])")
  print ("c.debug = '-g' # or other options for debugging purposes")
  print ("c.release = '-O3 -funroll-loops' # or other options for optimization")
  print ("c.compile = 'gcc -c -Wall ' # or other program for compiling")
  quit('Please add fortran compiler to config file')
  print ("cxx.ar = 'ar -r' # or other command for making library archives")

# Build individual libraries
try:
  database += ['blas', blas]
except:
  print( 'Add blas to config.py')
  sys.exit(1)

try:
  database += ['lapack', lapack, ['blas'] ]
except:
  print ('Add blas to config.py')
  sys.exit(1)

try:
  print ('ARPACK library: ', arpack.library)
  database += ['arpack', arpack, ['lapack']]
except:
  # Build ARPACK
  print ('ARPACK not found: building ARPACK')
  arpack = glas_buildsystem.library()
  arpack.library = '-L'+home+'/ARPACK/ -larpack' + ' ' + fortran.library
  database += ['arpack', arpack, ['lapack']]
  os.chdir('ARPACK')
  glas_test.write_make_dir( [], library_file )
  glas_test.make_dir(library_file, [])
  os.chdir('..')

#try:
#  print 'MCBSP include path: ', mcbsp.include_path
#  print 'MCBSP library: ', mcbsp.library
#  database += ['mcbsp', mcbsp]
#except:
#  # Build MCBSP
#  print 'MCBSP not found: building MCBSP'
#  if not 'cc' in locals():
#    print 'Help: I need the C compiler. Please add the following to config.py:'
#    print "cc = glas_buildsystem.new_compiler('C', 'c')"
#    print "cc.debug = '-g' # or other options for debugging purposes"
#    print "cc.release = '-O3 -funroll-loops -fvectorize' # or other options for optimization"
#    print "cc.compile = 'gcc -c -Wall ' # or program for compiling"
#    quit('Please add C compiler to config file')
#  mcbsp = glas_buildsystem.library()
#  mcbsp.include_path = [home+'/mcbsp']
#  mcbsp.library = '-L'+home+'/mcbsp/ -lmcbsp'
#  mcbsp.defines = '-DMCBSP_USE_MUTEXES -DMCBSP_COMPATIBILITY_MODE'
#  database += ['mcbsp', mcbsp]
#  os.chdir('mcbsp')
#  glas_test.write_make_dir( [], library_file )
#  glas_test.make_dir(library_file)
#  os.chdir('..')

# Boost
try:
  print ('Boost: ')
  print ('Boost include path: ') + boost.include_path
  database += ['boost', boost]
except:
  boost = glas_buildsystem.library()
  boost.include_path = [home+'/boost_1_59_0']
  boost.library = '-L'+home+'/boost_1_59_0/libs -lboost_1_59_0'
  database += ['boost',boost]

boost_numeric_bindings = glas_buildsystem.library()
boost_numeric_bindings.include_path = [home]
database += ['boost_numeric_bindings',boost_numeric_bindings, ['boost'] ]

arpack_bindings = glas_buildsystem.library()
arpack_bindings.include_path = [home+'/arpack_bindings']
arpack_bindings.library = '-L'+home+'/arpack_bindings/lib -larpack_bindings'
database += ['arpack_bindings',arpack_bindings ,['arpack','boost_numeric_bindings']]
os.chdir('arpack_bindings')
glas_test.write_make_dir( [], library_file )
glas_test.make_dir(library_file, [])
os.chdir('..')

mumps = glas_buildsystem.library()
mumps_root = home+'/MUMPS_5.0.1'
mumps.include_path = [mumps_root+'/include',mumps_root+'/libseq',mumps_root+'/pord/include',mumps_root+'/src']
mumps.library = '-L'+mumps_root+'/libseq -L'+mumps_root+'/lib -lcmumps -ldmumps -lmumps_common -lzmumps -lpord -lmpiseq ' + fortran.library
database += ['mumps',mumps , ['lapack']]
os.chdir(mumps_root)
glas_test.write_make_dir( [], library_file )
glas_test.make_dir(library_file, [])
os.chdir('..')

os.chdir('boost_1_59_0/libs')
glas_test.write_make_dir( [], library_file )
glas_test.make_dir(library_file, [])
os.chdir(home)

# Glas2
glas2 = glas_buildsystem.library()
glas2_root = home + 'glas2'
glas2.include_path = [glas2_root]
glas2.defines = '-DGLAS_COMPLEX'
glas2.library = '-L'+glas2_root+'/libs/ -lglas2'
database += ['glas2', glas2, ['boost']]

os.chdir(glas2_root+'/libs')
glas_test.write_make_dir( [], library_file )
glas_test.make_dir(library_file, [])
os.chdir(home)

# Itsol
itsol = glas_buildsystem.library()
itsol_root = home + '/ITSOL_2'
itsol.include_path = [itsol_root+'/INC']
itsol.library = '-L'+itsol_root + ' -litsol'
database += ['itsol', itsol]

os.chdir(itsol_root)
glas_test.write_make_dir( [], library_file )
glas_test.make_dir(library_file, [])
os.chdir(home)

itsol_pp = glas_buildsystem.library()
itsol_pp_root = home + '/itsol++'
itsol_pp.include_path = [itsol_pp_root]
itsol_pp.library = '-L'+itsol_pp_root+'/lib/ -litsol_pp'
database += ['itsol_pp', itsol_pp, ['itsol','boost_numeric_bindings'] ]

os.chdir('itsol++/lib')
glas_test.write_make_dir( [], library_file )
glas_test.make_dir(library_file, [])
os.chdir('../..')

#libraries['glas3'] = glas_buildsystem.library()
#glas3_root = home
#libraries['glas3'].include_path = [glas2_root]
#libraries['glas3'].library = '-L'+glas3_root+'/glas3/libs/ -lglas3'
#libraries['glas3'] += libraries['boost']
