#! /usr/bin/python

import glas_test

try:
  glas_test.make_dir()
except NameError:
  print 'Make failed'
