build = 'external'

def generate_make():
  f = open('Make_glas.inc', 'r')
  fw = open('Makefile.inc','w')
  fw.write(f.read())
  f.close()

  fw.write('#  GLAS additions' + '\n')

  fortran = find_compiler('Fortran')
  fw.write('FC = ' + fortran.compile + '\n' )
  fw.write('OPTF = ' + fortran.release + ' ' + fortran.flags + '\n' )

  cc = find_compiler('C')
  fw.write('CC = ' + cc.compile + '\n' )
  fw.write('OPTC = ' + cc.release + ' ' + cc.flags + '\n' )
  fw.write('CDEFS   = -D' + fortran.c_name_ext +'\n')
  fw.write('#  End GLAS additions' + '\n')
  fw.close()

do_make = 'make alllib'
