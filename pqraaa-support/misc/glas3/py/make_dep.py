#! /usr/bin/python

import os
import dep_tools
import sys

if len(sys.argv)>1:
  os.chdir( sys.argv[1] )
glas_root = os.getcwd()

make_f = open( os.path.join( os.getenv('HOME'), '.glas_build_dependencies.db' ), 'w' )
make_f.write('GLAS_ROOT = ' + glas_root + '\n\n' )

dict = {}
dep_tools.make_deps( dict, make_f, glas_root )
