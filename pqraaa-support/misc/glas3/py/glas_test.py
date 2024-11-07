import os
import sys
import dep_tools

def handle_dir(file,accum_dir):
  pwd = os.getcwd()
 # file.write('<tr><td colspan="3">' + pwd + '</td></tr>\n' )
  lst = os.listdir(pwd)
  lst.sort()

  # First do all executable files
  for f in lst:
    fn = os.path.split(f)[1] ;
    ind = fn.find('.cpp') ;
    if ind>0:
      exename = fn[0:ind]
      file.write('<tr><td>' + accum_dir + '</td><td><span style="color: blue"><a href="' + os.path.join(accum_dir,exename+'_log.txt') + '">' + exename + '</a></span></td>')
      if os.path.isfile( os.path.join('.',exename ) ):
        ret = os.system('./'+exename +' 2>&1 > ' + exename+'_log.txt')
        ret = ret / 256
        if ret==0:
          file.write('<td><span style="color: green">OK</span>' + '</td></tr>\n')
          print accum_dir + exename + ' : \033[32mOK\033[00m'
        else:
          file.write('<td><span style="color: red">' + str(ret) + '</span>' + '</td></tr>\n')
          print accum_dir + exename + ' : \033[31mKO\033[00m'
      else:
        print accum_dir + exename + ' : \033[33mNA\033[00m'
        file.write('<td><span style="color: orange">NA</span>' + '</td></tr>\n')

      file.write('<tr><td>' + accum_dir + '</td><td><span style="color: blue"><a href="' + os.path.join(accum_dir,exename+'_d_log.txt') + '">' + exename + '_d</a></span></td>')
      if os.path.isfile( os.path.join('.',exename+'_d' ) ):
        ret = os.system('./'+exename +'_d 2>&1 >' + exename+'_d_log.txt')
        ret = ret / 256
        if ret==0:
          file.write('<td><span style="color: green">OK</span>' + '</td></tr>\n')
          print accum_dir + exename + '_d : \033[32mOK\033[00m'
        else:
          file.write('<td><span style="color: red">' + str(ret) + '</span>' + '</td></tr>\n')
          print accum_dir + exename + '_d : \033[31mKO\033[00m'
#file.write('<tr><td>' + accum_dir + '</td><td><span style="color: blue">' + exename + '</span></td>')
      else:
        print accum_dir + exename + '_d : \033[33mNA\033[00m'
        file.write('<td><span style="color: orange">NA</span>' + '</td></tr>\n')


  # Then do all directories
  for f in lst:
    if os.path.isdir(f):
      if os.path.split(f)[1][0]!='.':
        os.chdir(f)
        handle_dir(file, accum_dir + os.path.split(f)[-1] + '/')
        os.chdir(pwd)

def clean_up_dir():
  pwd = os.getcwd()
  lst = os.listdir(pwd)
  lst.sort()

  # First do all executable files
  if os.path.exists( os.path.join('.', 'glas_build.py' ) ):
    for f in lst:
      fn = os.path.split(f)[1] ;
      ind = fn.find('.o') ;
      if ind>0:
        os.remove( f )
      ind = fn.find('.cpp') ;
      if ind>0:
        exename = fn[0:ind]
        if os.path.exists( exename ):
          os.remove( os.path.join('.', exename ) )
        if os.path.exists( exename+'_d' ):
          os.remove( os.path.join('.', exename+'_d' ) )
      ind = fn.find('_time') ;
      if ind>0:
        os.remove( os.path.join('.',fn ) )
      ind = fn.find('_log.txt')
      if ind>0:
        os.remove( os.path.join('.',fn ) )

  # Then do all directories
  for f in lst:
    if os.path.isdir(f):
      if os.path.split(f)[1][0]!='.':
        os.chdir(f)
        clean_up_dir()
        os.chdir(pwd)

def make_dir():
  #sys.path = sys.path + os.environ['HOME']
  #import make.config
  pwd = os.getcwd()
  lst = os.listdir(pwd)
  lst.sort()

  # First do all directories
  for f in lst:
    if os.path.isdir(f):
      if os.path.split(f)[1][0]!='.':
        os.chdir(f)
        make_dir()
        os.chdir(pwd)

  if not os.path.exists( os.path.join( os.getenv('HOME'), 'glas_build_init.py') ):
    return

  optf = open( os.path.join( os.getenv('HOME'), 'glas_build_init.py'), 'r' )
  for l in optf:
    exec( l )

  if os.path.exists( os.path.join('.', 'glas_build.py' ) ):
    # Read glas_build.py options file and execute line by line
    dependencies = glas_buildsystem.library()
    dependencies += glas3
    type = 'none'
    build = 'exe'
    variant = ['debug']
    try:
      optf = open( 'glas_build.py', 'r' )
      for l in optf:
        exec( l )
    except:
      print pwd + ' : An error ocurred: skip this directory'
      return
    print 'Compile directory "' + pwd + '"'
    if type!='none':
      build = type
      if type=='test':
        build = 'exe'
        variant = ['debug']
      elif type=='perf_test':
        build = 'exe'
        variant = ['release']


#if os.path.isfile( os.path.join('.', 'makefile' ) ):
    # First build all executable files
    if len(sys.argv)>1:
       for exename in sys.argv[1:]:
          if build=='lib' or build=='ar_lib' or build=='dyn_lib':
            print 'Compile obj ' + exename
            for v in variant:
              if v=='release':
                ret = os.system('make '+exename + '.o 2>&1 > ' + exename+'_log.txt' )
              elif v=='debug':
                ret = os.system('make '+exename + '_d.o 2>&1 > ' + exename+'_d_log.txt' )
          else:
            print 'Compile and link ' + exename
            for v in variant:
              if v=='release':
                ret = os.system('make '+exename + ' 2>&1 | tee ' + exename+'_log.txt' )
              elif v=='debug':
                ret = os.system('make '+exename + '_d 2>&1 | tee ' + exename+'_d_log.txt' )
    else:
      for f in lst:
        fn = os.path.split(f)[1] ;
        ind = fn.find('.cpp') ;
        exename = fn[0:ind]
        if ind>0:
          if os.path.exists( os.path.join( ".", exename ) ):
            os.remove( os.path.join( ".", exename ) )
          if os.path.exists( os.path.join( ".", exename + "_d" ) ):
            os.remove( os.path.join( ".", exename + "_d" ) )
          if build=='lib' or build=='ar_lib' or build=='dyn_lib':
            print 'Compile obj ' + exename
            for v in variant:
              if v=='release':
                ret = os.system('make '+exename + '.o 2>&1 > ' + exename+'_log.txt' )
              elif v=='debug':
                ret = os.system('make '+exename + '_d.o 2>&1 > ' + exename+'_d_log.txt' )
          else:
            print 'Compile and link ' + exename
            for v in variant:
              if v=='release':
                ret = os.system('make '+exename + ' 2>&1 | tee ' + exename+'_log.txt' )
              elif v=='debug':
                ret = os.system('make '+exename + '_d 2>&1 | tee ' + exename+'_d_log.txt' )
      if build=='ar_lib' or build=='dyn_lib':
        print 'Link ' + os.path.split(pwd)[1]
        ret = os.system('make 2>&1 > ar_log.txt')

class obj_list_type:
  def __init__( self ):
    self.debug = ''
    self.release = ''

  def add( self, that ):
    self.debug += ' ' + that.debug
    self.release += ' ' + that.release

def gen_make_ar_lib():
  pwd = os.getcwd()
  lst = os.listdir(pwd)
  lst.sort()

  obj_list = obj_list_type()

  # do all subdirectories
  for f in lst:
    if os.path.isdir(f):
      if os.path.split(f)[1][0]!='.':
        os.chdir(f)
        obj_list.add( gen_make_ar_lib() )
        os.chdir(pwd)

  optf = open( os.path.join( os.getenv('HOME'), 'glas_build_init.py'), 'r' )
  for l in optf:
    exec( l )

  if os.path.exists( os.path.join('.', 'glas_build.py' ) ):
    # Read glas_build.py options file and execute line by line
    dependencies = glas_buildsystem.library()
    dependencies += glas3
    try:
      optf = open( 'glas_build.py', 'r' )
      for l in optf:
        exec( l )
    except:
      print pwd + ' : An error ocurred: skip this directory'
      return obj_list
    else:
      for f in lst:
        fn = os.path.split(f)[1] ;
        ind = fn.find('.cpp') ;
        if ind>0:
          exename = fn[0:ind]
          obj_list.debug += ' ' + os.path.join( pwd, exename + '_d.o' )
          obj_list.release += ' ' + os.path.join( pwd, exename + '.o' )

  return obj_list

def write_make_dir( argv ):
  pwd = os.getcwd()
  lst = os.listdir(pwd)
  lst.sort()

  dict = {}

  if not os.path.exists( os.path.join( os.getenv('HOME'), 'glas_build_init.py') ):
    return

  optf = open( os.path.join( os.getenv('HOME'), 'glas_build_init.py'), 'r' )
  for l in optf:
    exec( l )

  # check if options are available
  if os.path.exists( os.path.join('.', 'glas_build.py' ) ):

    dependencies = glas_buildsystem.library()
    dependencies += glas3
    type = 'none'
    build = 'exe'
    variant = ['debug']

    # Read glas_build.py options file and execute line by line
    try:
      optf = open( 'glas_build.py', 'r' )
      for l in optf:
        exec( l )
    except:
      print pwd + ' : An error occurred: skip this directory'
    else:
      include_command = ''
      for i in glas_buildsystem.include_path(dependencies):
        include_command += ' -I' + i
      if type!='none':
        build = type
        if type=='test':
          build = 'exe'
          variant = ['debug']
        elif type=='perf_test':
          build = 'exe'
          variant = ['release']

      glas_roots = glas_buildsystem.modif_include_path(dependencies)

      cxxflags = glas_buildsystem.cxxflags( dependencies )
      library = glas_buildsystem.link( dependencies )
      defines = glas_buildsystem.defines( dependencies )
      ldflags = glas_buildsystem.ldflags( dependencies )

      print 'Generate makefile for directory "' + pwd + '" with build "' + build + '" and variants "' + str(variant) + '"'

      fm = open('makefile', 'w')
#fm.write('include ' + os.path.join( dep_tools.glas_root(), 'make_deps.db\n\n' ) )
      fm.write('CXXFLAGS = ' + cxx.cxxflags + ' ' + defines + ' ' + cxxflags + ' ' + include_command + '\n')
      fm.write('LDFLAGS = ' + cxx.ldflags + ' ' + ldflags + '\n')
      fm.write('LIBS = ' + library + ' ' + cxx.library + '\n')

      if build=='dyn_lib':
        libname_d = libname + '_d.so'
        libname = libname + '.so'
        obj_list = gen_make_ar_lib()

        fm.write('\n')
        fm.write('all:\n')
        for v in variant:
          if v=='release':
            fm.write('	@' + cxx.link + ' ' + libname + ' $(LDFLAGS) ' + cxx.shared_ldflags + ' ' + cxx.release + ' ' + obj_list.release + ' $(LIBS)\n' )
          elif v=='debug':
            fm.write('	@' + cxx.link + ' ' + libname_d + ' $(LDFLAGS) ' + cxx.shared_ldflags + ' ' + cxx.debug + ' ' + obj_list.debug + ' $(LIBS)\n' )
      elif build=='ar_lib':
        libname_d = libname + '_d.a'
        libname = libname + '.a'
        obj_list = gen_make_ar_lib()

        fm.write('\n')
        fm.write('all:\n')
        for v in variant:
          if v=='release':
            fm.write('	@' + cxx.ar + ' ' + libname + ' ' + obj_list.release + '\n' )
          elif v=='debug':
            fm.write('	@' + cxx.ar + ' ' + libname_d + ' ' + obj_list.debug + '\n' )
  
    #if os.path.isfile( os.path.join('.', 'makefile' ) ):
        # First build all executable files
      for f in lst:
          fn = os.path.split(f)[1] ;
          ind = fn.find('.cpp') ;
          if ind>0:
            exename = fn[0:ind]
            if fn==exename + '.cpp':
              fm.write('\n')
              deps_var = exename.upper() + '_DEPS'
              fm.write( deps_var + ' = ' )
              deps = dep_tools.make_file_deps( [], dict, [], glas_roots, f )
              for i in range(0,len(deps)):
                 fm.write( '\\\n	' + deps[i] )
              fm.write('\n\n')
              for v in variant:
                if v=='release':
                  fm.write(exename + '.o : ' + exename + '.cpp $(' + deps_var + ')\n' )
                  fm.write('	@' + cxx.compile + ' ' + cxx.release + ' $(CXXFLAGS) -o ' + exename + '.o ' + exename + '.cpp\n' )
                elif v=='debug':
                  fm.write(exename + '_d.o : ' + exename + '.cpp $(' + deps_var + ')\n')
                  fm.write('	@' + cxx.compile + ' ' + cxx.debug + ' $(CXXFLAGS) -o ' + exename + '_d.o ' + exename + '.cpp\n' )
                fm.write('\n')
              if build=='exe':
                for v in variant:
                  fm.write('\n')
                  if v=='release':
                    fm.write(exename + ': ' + exename + '.o\n')
                    fm.write('	@' + cxx.link + ' ' + exename + ' $(LDFLAGS) ' + exename + '.o $(LIBS)\n' )
                  elif v=='debug':
                    fm.write(exename + '_d: ' + exename + '_d.o\n')
                    fm.write('	@' + cxx.link + ' ' + exename + '_d $(LDFLAGS) ' + exename + '_d.o $(LIBS)\n' )

      fm.close()

  # Then do all directories
  for f in lst:
    if os.path.isdir(f):
      if os.path.split(f)[1][0]!='.':
        os.chdir(f)
        write_make_dir( argv )
        os.chdir(pwd)
