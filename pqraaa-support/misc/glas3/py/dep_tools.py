import os
import re

def glas_root():
  if not os.path.exists( os.path.join( os.getenv('HOME'), 'glas_build_init.py') ):
    return

  optf = open( os.path.join( os.getenv('HOME'), 'glas_build_init.py'), 'r' )
  for l in optf:
    exec( l )
  return glas.include_path.split()[0][2:]

def dep_transform( x ):
  x = x.replace('/','_')
  x = x.replace('.','_') ;
  x = x.upper() ;
  return x

def is_in( list, item ):
  for i in list:
    if item==i:
       return True
  return False

def match_file_name( glas_root, file ):
  for i in glas_root:
    if os.path.exists( os.path.join( i, file ) ):
      return os.path.join( i, file )
  return ''

def make_file_deps( deps, dict, done, glas_root, file ):
  if not os.path.exists( file ) :
    return deps
  f = open(file,'r')
  for l in f:
    ls = l.split()
    if len(ls)>0 and ls[0]=='#include':
      include_file = ls[1]
      if ls[1][0]=='<':
#full_file_name = os.path.join( glas_root, ls[1][1:-1] )
        full_file_name = match_file_name( glas_root, ls[1][1:-1] )
        if full_file_name!='':
          do_addition = True
          if not is_in( deps, full_file_name ):
            deps.append( full_file_name )
            if not dict.has_key( full_file_name ) :
              if is_in( done, full_file_name ):
                print 'Circular dependence detected - skip'
                do_addition = False
              else:
                done.append(full_file_name)
                dict.update( { full_file_name : make_file_deps( [], dict, done, glas_root, full_file_name ) } ) ;
          if do_addition:
            for i in dict[full_file_name]:
              if not is_in( deps, i ):
                deps.append( i )
      elif ls[1][0]=='"':
        full_file_name = os.path.join( os.getcwd(), os.path.split(file)[0], ls[1][1:-1] )
        if not is_in( deps, full_file_name ):
          deps.append( full_file_name ) 
          for i in make_file_deps( [], dict, done, glas_root, ls[1][1:-1] ):
            if not is_in( deps, i ):
               deps.append( i )
  return deps


def make_deps( dict, f_make, root_dir ):
  # Make dependencies in the files of the current directory.
  pwd = os.getcwd()
  print 'Make dependencies in ' + pwd

  # Loop over files in the directory and add dependencies to the dictionary d.
  lst = os.listdir(pwd)
  lst.sort()
  for f in lst:
    fn = os.path.split(f)[1] ;
    ind = fn.find('.hpp') ;
    if ind>0 and fn[-1]=='p':
      f_make.write( '\n' + dep_transform( os.path.join( root_dir, f ) ) + ' = ' + os.path.join( pwd, f ) + '\\\n' )
      deps = make_file_deps( [], dict, false, root_dir, f )
      for names in deps:
        f_make.write( name + '\\\n' )
      dict.update( { os.path.join( root_dir, f ) : deps } )

  for f in lst:
    if os.path.isdir(f):
      if os.path.split(f)[1][0]!='.':
        os.chdir(f)
        dict = make_deps( dict, f_make, os.path.join( root_dir, f ) )
        os.chdir(pwd)
  return dict
