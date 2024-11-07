#! /usr/bin/python

import sys
import os
import datetime

def process_test( file, f ):
  if os.path.isfile( f ):
    fp = open( f )
    line = fp.readline().split() ;
    if len(line)>1:
      file.write('<td bgcolor="lightblue">' + line[1] + '</td>' ) ;
    else:
      file.write('<td bgcolor="lightblue">NA</td>' ) ;
    line = fp.readline().split() ;
    if len(line)>1:
      file.write('<td bgcolor="lightgreen">' + line[1] + '</td>' ) ;
    else:
      file.write('<td bgcolor="lightblue">NA</td>' ) ;
    line = fp.readline().split() ;
    if len(line)>1:
      file.write('<td bgcolor="yellow">' + line[1] + '</td>' ) ;
    else:
      file.write('<td bgcolor="lightblue">NA</td>' ) ;
  else:
    file.write('<td></td><td></td><td></td>')


def process_dir( name, case_list ):
  lst = os.listdir('.')
  lst.sort()
  old_dir = os.getcwd()
  for f in lst:
    if os.path.isdir(f):
      os.chdir( f )
      process_dir( name + ' ' + str(f), case_list )
      os.chdir( old_dir )

  do_it = 0
  for case in case_list:
    for case in case_list:
      if os.path.exists( os.path.join('.',case) ):
        do_it = 1

  if do_it:
    file.write('<tr><td>' + name + '</td>')
    for case in case_list:
      if os.path.exists( os.path.join('.',case) ):
        time_file = case + '_time'
        if ( os.path.exists( os.path.join('.', time_file ) ) ):
          process_test( file, os.path.join('.', time_file ) )
        else:
          print 'There is no "' + time_file + '"'
    file.write('</tr>\n')
  os.chdir( old_dir )

case_list = ['glas', 'blas', 'ublas']
if len(sys.argv)>1:
  case_list = [];
  for i in range(1,len(sys.argv)):
    case_list.append( sys.argv[i] )
print 'Cases : ' + str( case_list )

file = open('report.html', 'w') ;
file.write('<body>\n')
file.write('<title>Glas: Performance Test report</title>\n')
file.write('<h1>Report on GLAS performance tests (' + str(datetime.date.today()) + ')</h1>\n')
file.write('<table border="1"><tr><th></th>')
for case in case_list:
  file.write('<th colspan="3">' + case + '</th>')
file.write('</tr>\n')
file.write('<tr><th>Test case:</th>')
for case in case_list:
  file.write('<td bgcolor="lightblue">100</td><td bgcolor="lightgreen">10000</td><td bgcolor="yellow">1000000</td>')
file.write('</tr>\n')

process_dir('', case_list)

file.write('</table>\n')
file.write('Timings in micro seconds\n')
file.write('</body>\n')
