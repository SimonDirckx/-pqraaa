#! /usr/bin/python

import sys
import os

def handle_dir(indent):
  wd = os.getcwd()

  if os.path.exists( os.path.join('.', 'glas' ) ):
     print indent + ' glas'
     ret = os.system( './glas > glas_log.txt' )
  if os.path.exists( os.path.join('.', 'blas' ) ):
     print indent + ' blas'
     ret = os.system( './blas > blas_log.txt' )
  if os.path.exists( os.path.join('.', 'ublas' ) ):
     print indent + ' ublas'
     ret = os.system( './ublas > ublas_log.txt' )

  lst = os.listdir('.')
  lst.sort()
  for f in lst:
    if os.path.isdir(f):
      if os.path.split(f)[1][0]!='.':
        print indent + 'Running : ' + f + ': '
        os.chdir( f )
        handle_dir( indent + '  ' )
        os.chdir( wd )

handle_dir('')
