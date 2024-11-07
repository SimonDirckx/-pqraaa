import glas_buildsystem

class gcc_compiler( glas_buildsystem.compiler ):
  def __init__(self):
    self.bin = 'g++'
    self.debug = '-g'
    self.release = '-O3 -DNDEBUG'

  def __str__(self):
    return 'gcc'
