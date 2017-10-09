################################################################################
#                                                                              #
#  CONFIGURATION AND COMPILATION ROUTINE                                       #
#                                                                              #
################################################################################

import sys
import os
import util
from subprocess import call
import subprocess
      
NOPARAM = False
PARAM_NAME = 'param_template.dat'

CPARMS = {}
RPARMS = {}

# CONTAINER CLASS FOR RUNTIME PARAMETER
class rparm:
  datatype = None
  value = None

# SET COMPILE-TIME PARAMETER
def set_cparm(name, value):
  CPARMS[name] = value

# SET RUNTIME PARAMETER. DO NOT OVERWRITE DEFAULT VALUES, AS THIS IS CALLED BY
# PROBLEM FILE BEFORE CORE ROUTINE
def set_rparm(name, datatype, default=None):
  if datatype == 'int':
    datatype = 'integer'
  
  if datatype != 'integer' and datatype != 'double' and datatype != 'string':
    util.warn('DATATYPE ' + datatype + ' NOT SUPPORTED')
    print 'CHOICES: integer double string'
    sys.exit()
  
  obj = rparm()
  obj.datatype = datatype

  if name in RPARMS:
    obj.value = RPARMS[name].value
  else:
    obj.value = default
  RPARMS[name] = obj

def write_rparm(pf, name):
  if name not in RPARMS:
    util.warn('RUNTIME PARAMETER ' + name + ' NOT SET')
    sys.exit()

  datatype = RPARMS[name].datatype
  default = RPARMS[name].value
  if datatype == 'integer':
    if default == None:
      pf.write('[int] ' + name + ' = \n')
    else:
      pf.write('[int] ' + name + ' = %d\n' % default)
  elif datatype == 'double':
    if default == None:
      pf.write('[dbl] ' + name + ' = \n')
    else:
      pf.write('[dbl] ' + name + ' = %e\n' % default)
  elif datatype == 'string':
    if default == None:
      pf.write('[str] ' + name + ' = \n')
    else:
      pf.write('[str] ' + name + ' = %s\n' % default)
  else:
    print name
    util.warn("DATATYPE " + str(datatype) + " NOT RECOGNIZED")
    sys.exit()

  del RPARMS[name]

def is_user_debug():
  for n in xrange(len(sys.argv)):
    if (sys.argv[n] == '-debug'):
      return 1
  return 0

def set_dirs(PATHS):
  MOVEEXEC = True
  PATHS['BUILD'] = 'SRC'
  PATHS['PROB'] = os.getcwd()
  for n in xrange(len(sys.argv)):
    if sys.argv[n] == '-dir' and n < len(sys.argv) - 1:
      MOVEEXEC = False
      PATHS['BUILD'] = util.sanitize_path(sys.argv[n+1])
  PATHS['SRC'] = os.path.join(PATHS['BUILD'], 'source')

  for key in PATHS:
    if PATHS[key][0] != '/':
      PATHS[key] = os.path.join(PATHS['PROB'], PATHS[key], '')
    util.make_dir(PATHS[key])
  return MOVEEXEC

def build(PROBLEM, PATHS):
  # PROCESS USER INPUT
  DEBUG = is_user_debug()
  MOVEEXEC = set_dirs(PATHS)

  # PRINT TO TERMINAL AND LOGFILE
  LOGFILE = os.path.join(PATHS['BUILD'], 'log_build')
  util.log_output(sys, LOGFILE)

  for n in xrange(len(sys.argv)):
    if sys.argv[n] == '-noparam':
      NOPARAM = 1

  print("")
  print("********************************************************************************")
  print("")
  print("                              BHLIGHT BUILD SCRIPT")
  print("")
  print("  OPTIONS:")
  print("    -debug               (compile files sequentially)")
  print("    -noparam             (do not create new parameter file)")
  print("    -dir /path/to/target (create target dir and build there)")
  print("")
  print("********************************************************************************")
  print("")

  # SEARCH FOR MACHINE
  machines = util.get_files(PATHS['MACHINE'], '*')
  for n in xrange(len(machines)):
    machines[n] = machines[n].split('/')[-1].replace('.py', '')
    machine = __import__(machines[n])
    if machine.matches_host() == True:
      break
    del machine

  try: machine
  except NameError: util.warn("HOST " + os.uname()[1] + " UNKNOWN"); sys.exit()

  host = machine.get_options()

  C_FLAGS = '-std=c99 -mcmodel=medium ' + host['COMPILER_FLAGS']

  # MATH AND DYNAMIC LINKING
  LIB_FLAGS = '-lm -ldl'

  LIBRARIES = ''
  INCLUDES  = ''

  # GSL
  host['GSL_DIR'] = util.sanitize_path(host['GSL_DIR'])
  LIB_FLAGS += ' -lgsl -lgslcblas'
  LIBRARIES += '-L' + host['GSL_DIR'] + 'lib/'
  INCLUDES  += '-I' + host['GSL_DIR'] + 'include/'

  # MPI
  if 'MPI_DIR' in host:
    host['MPI_DIR'] = util.sanitize_path(host['MPI_DIR'])
    LIBRARIES += ' -L' + host['MPI_DIR'] + 'lib/'
    INCLUDES  += ' -I' + host['MPI_DIR'] + 'include/'

  # HDF5
  if 'HDF5_DIR' in host:
    host['HDF5_DIR'] = util.sanitize_path(host['HDF5_DIR'])
    LIB_FLAGS += ' -lhdf5_hl -lhdf5'
    LIBRARIES += ' -L' + host['HDF5_DIR'] + 'lib/'
    INCLUDES  += ' -I' + host['HDF5_DIR'] + 'include/'

  print("  CONFIGURATION\n")
  def print_config(key, var):
    print("    " + util.color.BOLD + "{:<15}".format(key) + util.color.NORMAL +
          str(var))

  print_config("MACHINE",    host['NAME'])
  print_config("PROBLEM",    PROBLEM)
  print_config("COMPILER",   host['COMPILER'])
  print_config("GSL_DIR",    host['GSL_DIR'])
  if 'MPI_DIR' in host:
    print_config("MPI_DIR",  host['MPI_DIR'])
  if 'HDF5_DIR' in host:
    print_config("HDF5_DIR", host['HDF5_DIR'])
  print_config("EXECUTABLE", host['EXECUTABLE'])
  print_config("C_FLAGS",    C_FLAGS)
  print_config("LIB_FLAGS",  LIB_FLAGS)
  print_config("LIBRARIES",  LIBRARIES)
  print_config("INCLUDES",   INCLUDES)
  print_config("OPENMP", CPARMS['OPENMP'])

  print("\n  COMPILE-TIME PARAMETERS\n")
  print_config("N1TOT", CPARMS['N1TOT'])
  print_config("N2TOT", CPARMS['N2TOT'])
  print_config("N3TOT", CPARMS['N3TOT'])
  print_config("N1CPU", CPARMS['N1CPU'])
  print_config("N2CPU", CPARMS['N2CPU'])
  print_config("N3CPU", CPARMS['N3CPU'])
  print_config("METRIC", CPARMS['METRIC'])
  print_config("RECONSTRUCTION", CPARMS['RECONSTRUCTION'])
  if util.parm_is_active(CPARMS, 'RADIATION'):
    print_config("RADIATION", CPARMS['RADIATION'])
  if util.parm_is_active(CPARMS, 'ELECTRONS'):
    print_config("ELECTRONS", CPARMS['ELECTRONS'])

  # Set core runtime parameters
  set_rparm('tf', 'double')
  set_rparm('dt', 'double')
  if CPARMS['METRIC'] == 'MINKOWSKI':
    set_rparm('x1Min', 'double', default = 0.)
    set_rparm('x1Max', 'double', default = 1.)
    set_rparm('x2Min', 'double', default = 0.)
    set_rparm('x2Max', 'double', default = 1.)
    set_rparm('x3Min', 'double', default = 0.)
    set_rparm('x3Max', 'double', default = 1.)
  if CPARMS['METRIC'] == 'MKS':
    set_rparm('a', 'double', default = 0.5)
    set_rparm('hslope', 'double', default = 0.3)
    set_rparm('poly_xt', 'double', default = 0.82)
    set_rparm('poly_alpha', 'double', default = 14.)
    set_rparm('mks_smooth', 'double', default = 0.5)
    set_rparm('Rout', 'double', default = 40.)
    #if util.parm_is_active(CPARMS, 'RADIATION'):
    #  set_rparm('Rout_rad', 'double')

  if (util.parm_is_active(CPARMS, 'RADIATION') or                                
      util.parm_is_active(CPARMS, 'COULOMB')):
    if CPARMS['METRIC'] == 'MINKOWSKI':
      set_rparm('L_unit', 'double')
      set_rparm('M_unit', 'double')
    if CPARMS['METRIC'] == 'MKS':
      set_rparm('M_unit', 'double')
      set_rparm('mbh', 'double', default = 1.989e34)
  
  set_rparm('gam', 'double', default = 5./3.)
  set_rparm('cour', 'double', default = 0.4)
  
  if util.parm_is_active(CPARMS, 'ELECTRONS'):
    set_rparm('game', 'double', default = 4./3.)
    set_rparm('gamp', 'double', default = 5./3.)
    set_rparm('fel0', 'double', default = 0.01)
    set_rparm('tptemin', 'double', default = 1.e-3)
    set_rparm('tptemax', 'double', default = 1.e3)
  
  if util.parm_is_active(CPARMS, 'RADIATION'):
    if not util.parm_is_active(CPARMS, 'ELECTRONS'):
      set_rparm('tp_over_te', 'double', default = 1.)
    set_rparm('nph_per_proc', 'double', default = 1.e5)
    set_rparm('numin', 'double', default = 1.e8)
    set_rparm('numax', 'double', default = 1.e20)
    set_rparm('tune_emiss', 'double', default = 1.)
    set_rparm('tune_scatt', 'double', default = 1.)
    set_rparm('thetae_max', 'double', default = 1.e3)
    set_rparm('sigma_max', 'double', default = 1.)
    set_rparm('kdotk_tol', 'double', default = 1.e-6)
    set_rparm('Nph_to_track', 'double', default = 1.e4);
  
  set_rparm('DTd', 'double', default = 0.5)
  set_rparm('DTl', 'double', default = 0.1)
  set_rparm('DTr', 'integer', default = 1024)
  set_rparm('DTp', 'integer', default = 100)
  set_rparm('DTf', 'integer', default = 10)
  set_rparm('outputdir', 'string', default = './')

  # GET ALL SOURCE FILES
  SRC_CORE = util.get_files(PATHS['CORE'], '*.c')
  INC_CORE = util.get_files(PATHS['CORE'], '*.h')
  SRC_PROB = util.get_files(PATHS['PROB'], '*.c')
  INC_PROB = util.get_files(PATHS['PROB'], '*.h')
  SRC = []
  INC = []

  # COPY SOURCE FILES TO BUILD_DIR
  for n in xrange(len(SRC_CORE)):
    call(['cp', SRC_CORE[n], PATHS['SRC'] + SRC_CORE[n].rsplit('/',1)[1]])
  for n in xrange(len(INC_CORE)):
    call(['cp', INC_CORE[n], PATHS['SRC'] + INC_CORE[n].rsplit('/',1)[1]])
  for n in xrange(len(SRC_PROB)):
    call(['cp', SRC_PROB[n], PATHS['SRC'] + SRC_PROB[n].rsplit('/',1)[1]])
  for n in xrange(len(INC_PROB)):
    call(['cp', INC_PROB[n], PATHS['SRC'] + INC_PROB[n].rsplit('/',1)[1]])

  # WRITE PARAMETERS FILE
  #pf = open(PATHS['SRC'] + 'params.h', 'w')
  #for KEY in CPARMS:
  #  if isinstance(CPARMS[KEY], str):
  #    pf.write("#define " + KEY + " (" + CPARMS[KEY] + ")\n")
  #  else:
  #    pf.write("#define " + KEY + " (%g)\n" % CPARMS[KEY])
  #pf.close()

  # GET SINGLE LISTS OF ALL SOURCE, OBJECT, AND HEADER FILES
  SRC_ALL = util.get_files(PATHS['SRC'], '*.c')
  INC_ALL = util.get_files(PATHS['SRC'], '*.h')
  SRC = ''
  OBJ = ''
  INC = ''
  for n in xrange(len(SRC_ALL)):
    SRC += '%s ' % SRC_ALL[n]
    OBJ += '%s.o ' % SRC_ALL[n][:-2]
  for n in xrange(len(INC_ALL)):
    INC += '%s ' % INC_ALL[n]

  # WRITE MAKEFILE
  os.chdir(PATHS['SRC'])
  mf = open('makefile', 'w')
  mf.write('CC = ' + host['COMPILER'] + '\n')
  mf.write('CCFLAGS = ' + C_FLAGS + ' ' + LIBRARIES + ' ' + INCLUDES + '\n')
  mf.write('LIB_FLAGS = ' + LIB_FLAGS + '\n')
  mf.write('CC_COMPILE = $(CC) $(CCFLAGS) -c' + '\n')
  mf.write('CC_LOAD = $(CC) $(CCFLAGS)' + '\n')
  mf.write('.c.o:' + '\n')
  mf.write('\t$(CC_COMPILE) $*.c' + '\n')
  mf.write('EXE = bhlight' + '\n')
  mf.write('all: $(EXE)' + '\n')
  mf.write('SRC = ' + SRC + '\n')
  mf.write('OBJ = ' + OBJ + '\n')
  mf.write('INC = ' + INC + '\n')
  mf.write('$(OBJ) : $(INC) makefile' + '\n')
  mf.write('$(EXE): $(OBJ) $(INC) makefile' + '\n')
  mf.write('\t$(CC_LOAD) $(OBJ) $(LIB_FLAGS) -o $(EXE)')
  mf.close()

  print("\n  COMPILING SOURCE\n")

  ncomp = 0
  first_error = 1
  if DEBUG:
    popen = subprocess.Popen(['make'], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, universal_newlines=True)
  else:
    popen = subprocess.Popen(['make','-j','10'], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, universal_newlines=True)

  for stdout_line in iter(popen.stdout.readline, ""):
    if stdout_line.rstrip()[-2:] == '.c':
      print("    [" + util.color.BOLD + util.color.BLUE +
            "%2d%%" % (100.*float(ncomp)/len(SRC_ALL)) + util.color.NORMAL +
            "] " + util.color.BOLD +
            stdout_line.rsplit('-c',1)[1].rstrip().lstrip().split('/')[-1] +
            util.color.NORMAL)
      ncomp += 1
  for stderr_line in iter(popen.stderr.readline, ""):
    # THIS ALSO FAILS FOR WARNINGS!!!
    if first_error == 1:
      util.warn("COMPILER ERROR")
      first_error = 0
    print(stderr_line.rstrip())

  if first_error != 1:
    util.warn("COMPILATION FAILED")
    sys.exit()

  obj_files = util.get_files(PATHS['SRC'], '*.o')
  for f in obj_files:
    os.remove(f)
  os.rename(PATHS['SRC'] + 'bhlight', PATHS['BUILD'] + 'bhlight')

  print("\n  BUILD SUCCESSFUL")

  # CREATE RUNTIME PARAMETERS FILE
  PARAMFILE = PATHS['BUILD'] + PARAM_NAME
  pf = open(PARAMFILE, 'w')
  pf.write("### RUNTIME PARAMETERS ###\n")
  pf.write("\n# COORDINATES\n")
  write_rparm(pf, 'tf')
  write_rparm(pf, 'dt')
  
  if CPARMS['METRIC'] == 'MINKOWSKI':
    write_rparm(pf, 'x1Min')
    write_rparm(pf, 'x1Max')
    write_rparm(pf, 'x2Min')
    write_rparm(pf, 'x2Max')
    write_rparm(pf, 'x3Min')
    write_rparm(pf, 'x3Max')
  if CPARMS['METRIC'] == 'MKS':
    write_rparm(pf, 'Rout')
    if util.parm_is_active(CPARMS, 'RADIATION'):
      write_rparm(pf, 'Rout_rad')
  
  if (util.parm_is_active(CPARMS, 'RADIATION') or
      util.parm_is_active(CPARMS, 'COULOMB')):
    pf.write("\n# UNITS\n")
    if CPARMS['METRIC'] == 'MINKOWSKI':
      write_rparm(pf, 'L_unit')
      write_rparm(pf, 'M_unit')
    if CPARMS['METRIC'] == 'MKS':
      write_rparm(pf, 'mbh')
      write_rparm(pf, 'M_unit')
  
  pf.write("\n# FLUID\n")
  write_rparm(pf, 'gam')
  write_rparm(pf, 'cour')
  
  if util.parm_is_active(CPARMS, 'ELECTRONS'):
    pf.write("\n# ELECTRONS\n")
    write_rparm(pf, 'game')
    write_rparm(pf, 'gamp')
    write_rparm(pf, 'fel0')
    write_rparm(pf, 'tptemin')
    write_rparm(pf, 'tptemax')
  
  if util.parm_is_active(CPARMS, 'RADIATION'):
    pf.write("\n# RADIATION\n")
    if not util.parm_is_active(CPARMS, 'ELECTRONS'):
      write_rparm(pf, 'tp_over_te')
    write_rparm(pf, 'nph_per_proc')
    write_rparm(pf, 'numin')
    write_rparm(pf, 'numax')
    write_rparm(pf, 'tune_emiss')
    write_rparm(pf, 'tune_scatt')
    write_rparm(pf, 'thetae_max')
    write_rparm(pf, 'sigma_max')
    write_rparm(pf, 'kdotk_tol')
    write_rparm(pf, 'Nph_to_track')
  
  pf.write("\n# OUTPUT\n")
  write_rparm(pf, 'DTd')
  write_rparm(pf, 'DTl')
  write_rparm(pf, 'DTr')
  write_rparm(pf, 'DTp')
  write_rparm(pf, 'outputdir')
  if len(RPARMS.keys()) > 0:
    pf.write("\n# PROBLEM\n")
    prob_keys = RPARMS.keys()
    for key in prob_keys:
      write_rparm(pf, key)
  pf.close()
  print("\n  RUNTIME PARAMETER FILE CREATED")

  if MOVEEXEC:
    os.rename(PATHS['BUILD'] + 'bhlight',   PATHS['BUILD'] + '../bhlight')
    os.rename(PATHS['BUILD'] + PARAM_NAME, PATHS['BUILD'] + '../' + PARAM_NAME)

  print("")

  sys.exit()

