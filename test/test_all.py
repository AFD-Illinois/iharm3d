################################################################################ 
#                                                                              # 
# RUN ALL TESTS AND CHECK FOR ACCURACY                                         # 
#                                                                              # 
################################################################################

import os
import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../script/')
import util
import subprocess as sp
import numpy as np
import glob
import pickle
from scipy.interpolate import interp1d as interp
import time

# SPECIFY EMAIL OPTIONS
TO = ['bryan10@illinois.edu']
CC = []
FROM = 'afd.illinois.testing@gmail.com'
PASSWORD = 'whatpasswordshouldichoose'
SUBJECT = 'BHLIGHT TESTING REPORT'
LOGNAME = 'test_all.txt'

# REPORT TIME FOR TESTING

# GET ALL TEST SCRIPTS
TESTS = glob.glob('*.py')
TESTS.remove('test_all.py')

SEND_REPORT = False
for arg in sys.argv:
  if arg == '-email':
    SEND_REPORT = True

print("")                                                                      
print("********************************************************************************")
print("")                                                                      
print("                                AUTOMATED TESTING")                    
print("")                                                                      
print("********************************************************************************")

util.log_output(sys, LOGNAME)

DATE = time.strftime('%Y/%m/%d')
TIME = time.strftime('%H:%M:%S')
MACHINE = os.uname()[1]
popen = sp.Popen(['git', 'show', '-s', '--format=%H'], stdout=sp.PIPE,
                   universal_newlines=True)
for line in iter(popen.stdout.readline, ""):
  HASH = line.lstrip().rstrip()
popen = sp.Popen(['git', 'branch'], stdout=sp.PIPE, universal_newlines=True)
for line in iter(popen.stdout.readline, ""):
  if line[0] == '*': 
    BRANCH = line[2:].rstrip()
print '\n  DATE:    ' + DATE
print '  TIME:    ' + TIME
print '  MACHINE: ' + MACHINE
print '  BRANCH:  ' + BRANCH
print '  COMMIT:  ' + HASH + '\n'

# USE INTERPOLATION ON A (ANALYTIC SOLUTION) TO COMPARE TO B
def L1_norm(xa, ya, xb, yb):
  fa = interp(xa, ya)
  norm = 0.
  for n in xrange(len(xb)):
    norm += np.fabs(yb[n] - fa(xb[n]))/((yb[n] + fa(xb[n]))/2.)

  return (norm/n)

FAIL = False
for TEST in TESTS:
  if TEST != 'sod.py' and TEST != 'brem.py':
    continue
  print '  ' + util.color.BOLD + TEST + util.color.NORMAL
  popen = sp.Popen(['python', TEST, '-auto'], stdout=sp.PIPE, stderr=sp.PIPE,
                   universal_newlines=True)
  for line in iter(popen.stdout.readline, ""):
    if line.lstrip().rstrip() == 'BUILD SUCCESSFUL':
      print '    BUILD SUCCESSFUL'  
  print '    RUN FINISHED' 
  
  data = pickle.load(open('data.p', 'rb'))

  norm = L1_norm(data['SOL'][0], data['SOL'][1], 
                 data['CODE'][0], data['CODE'][1])

  print '    ERROR: %.2g %%' % (100*norm)
  if norm < 0.01:
    print util.color.BOLD + '    PASS' + util.color.NORMAL + '\n'
  else:
    print util.color.WARNING + '    FAIL' + util.color.NORMAL + '\n'
    FAIL = True

  sp.call(['rm', 'data.p'])
 
if not SEND_REPORT:
  sp.call(['rm', LOGNAME])
  sys.exit()

import smtplib
if FAIL:
  SUBJECT += ' - FAIL'
else:
  SUBJECT += ' - PASS'
MESSAGE = ''
MFILE = open(LOGNAME, 'rb')
for line in MFILE:
  MESSAGE += line
MFILE.close()
EMAIL = ('From: %s\r\n' % FROM
         + 'To: %s\r\n' % ','.join(TO)
         + 'CC: %s\r\n' % ','.join(CC)
         + 'Subject: %s\r\n' % SUBJECT
         + '\r\n'
         + MESSAGE)
ADDRS = TO + CC
srvr = smtplib.SMTP('smtp.gmail.com', 587)
srvr.ehlo()
srvr.starttls()
srvr.ehlo()
srvr.login(FROM, PASSWORD)
srvr.sendmail(FROM, ADDRS, EMAIL)
srvr.close()
sp.call(['rm', LOGNAME])

