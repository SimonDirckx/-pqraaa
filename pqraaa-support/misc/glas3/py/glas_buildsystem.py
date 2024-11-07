def isin( list, x ):
  for i in list:
    if i==x:
      return True
  return False

def merge( y, x ):
  for i in x:
    if not isin( y, i ):
      y.append( i )
  return y

class library:
  def __init__(self):
    self.include_path = []
    self.library = ''
    self.cxxflags = ''
    self.ldflags = ''
    self.defines = ''
    self.dependencies = []
    self.check_modif = False

  def __iadd__(self,that):
#    self.include_path = merge( self.include_path, that.include_path )
#    self.library += ' ' + that.library
#    self.cxxflags += ' ' + that.cxxflags
#    self.ldflags += ' ' + that.ldflags
#    self.defines += ' ' + that.ldflags
    if not isin( self.dependencies, that ):
      self.dependencies.append( that )
    self.dependencies = merge( self.dependencies, that.dependencies )
    return self

  def __str__(self):
    s = 'Library\n'
    s += ' include_path = ' + str(self.include_path) + '\n'
    s += ' library = ' + self.library+ '\n'
    s += ' cxxflags = ' + self.cxxflags+ '\n'
    s += ' ldflags = ' + self.ldflags+ '\n'
    s += ' defines = ' + self.defines
    s += ' check_modif = ' + str(self.check_modif)
    return s

def include_path( deps ):
  l = deps.include_path
  for i in deps.dependencies:
     l = merge( l, i.include_path )
  return l

def link( deps ):
  l = deps.library
  for i in deps.dependencies:
     l += ' ' + i.library
  return l

def cxxflags( deps ):
  l = deps.cxxflags
  for i in deps.dependencies:
     l += ' ' + i.cxxflags
  return l

def ldflags( deps ):
  l = deps.ldflags
  for i in deps.dependencies:
     l += ' ' + i.ldflags
  return l

def defines( deps ):
  l = deps.defines
  for i in deps.dependencies:
     l += ' ' + i.defines
  return l

def modif_include_path( deps ):
  if deps.check_modif:
    l = deps.include_path
  else:
    l = []
  for i in deps.dependencies:
     if i.check_modif:
       l = merge( l, i.include_path )
  return l

class compiler:
  def __init__(self):
    self.debug = '-g'
    self.release = '-O -DNDEBUG'
    self.compile = ''
    self.link = ''
    self.include_path = []
    self.include_dir_prefix = '-I'
    self.link_dir_prefix = '-L'
    self.cxxflags = ''
    self.ldflags = ''
    self.library = ''
    self.shared_cxxflags = ''
    self.shared_ldflags = ''
    self.ar = 'ar -rv'

  def __str__(self):
    s = 'Compiler\n'
    s += ' debug = ' + self.debug+ '\n'
    s += ' release = ' + self.release+ '\n'
    s += ' compile = ' + self.compile+ '\n'
    s += ' link = ' + self.link+ '\n'
    s += ' include_path = ' + self.include_path+ '\n'
    s += ' include_dir_prefix = ' + self.include_dir_prefix+ '\n'
    s += ' link_dir_prefix = ' + self.link_dir_prefix+ '\n'
    s += ' cxxflags = ' + self.cxxflags+ '\n'
    s += ' ldflags = ' + self.ldflags+ '\n'
    s += ' library = ' + self.library+ '\n'
    s += ' shared_cxxflags = ' + self.shared_cxxflags+ '\n'
    s += ' shared_ldflags = ' + self.shared_ldflags
    s += ' ar = ' + self.ar
    return s
