#! /usr/bin/python

import sys
import os
import datetime

import glas_test

file = open('report.html', 'w') ;
file.write('<body>\n')
file.write('<title>Glas: Test report</title>\n')
file.write('<h1>Report on GLAS regression tests (' + str(datetime.date.today()) + ')</h1>\n')
file.write('<table><tr><th>Directory</th><th>Test</th><th>Result</th></tr>\n')
glas_test.handle_dir(file,'./')
file.write('</table>\n')
file.write('</body>\n')
