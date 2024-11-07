#! /usr/bin/python

import sys
import os
import shutil

def copy_dirs( fr, to ):
  if not os.path.exists(to):
    os.mkdir(to)
  for l in os.listdir(fr):
    if l.count('.svn')==0:
      if os.path.isdir(fr + '/'+l):
        if not os.path.exists( to + '/' + l ):
           os.mkdir( to + '/' + l )
        copy_dirs( fr + '/' + l, to + '/' + l )
      else:
        print fr + '/'+l + ' ==> ' + to + '/' + l
        shutil.copyfile( fr +'/'+l, to + '/' + l )

if len(sys.argv)<=2:
  print 'Usage: copy_sources <from> <to>'
else:
  from_dir = sys.argv[1]
  to_dir = sys.argv[2]

  copy_dirs( from_dir, to_dir )
