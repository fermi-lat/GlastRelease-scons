#!/usr/bin/env python
import sys, string, os, re ;


#=================================================
# globals
#=================================================

original_dir = os.getcwd()
release_expr = re.compile('^v[0-9]+r[0-9]+')


#=================================================
# the function below establishes the path
# to the root directory of a given client package
#=================================================

def cmt_path(package) :

  # check if it is the original package
  if original_dir.rfind(os.path.join('',package,'')) != -1  :
    return os.path.join('..','..','cmt')
  
  # ask cmt
  packages_pipe = os.popen('cmt show packages')
  for line in packages_pipe :
    tokens = line.split()
    if tokens[0] == package :
      return os.path.join(tokens[2],tokens[0],tokens[1],'cmt')
  packages_pipe.close()

  # not found
  print 'VALIDATION ERROR: package',package,'NOT FOUND'
  sys.exit(1)
      

#=================================================
# build a package test application
#=================================================

def build_application_test(package) :

  os.chdir(cmt_path(package))
  
  build_command = 'cmt bro -local cmt config'
  if os.name == 'posix':
    build_command += ' ; cmt bro -local make'
    build_command += ' ; make test'
    if os.system(build_command) != 0 :
      print 'VALIDATION ERROR: test_'+package+'.exe BUILD FAILED'
      sys.exit(1)
# David: Windows nmake compilation fails for some reason,
# so one will need to compile interactively with MRvcmt
# before launching this validation script
#
#  if os.name == 'nt':
#    build_command += ' & cmt bro -local nmake /f nmake'
#    build_command += ' & nmake /f nmake test'
#    if os.system(build_command) != 0 :
#      print 'VALIDATION ERROR: test_'+package+'.exe BUILD FAILED'
#      sys.exit(1)

  os.chdir(original_dir)


#=================================================
# help output analysis
#=================================================

def key_values(filename) :

  # lines of interest
  regexp = re.compile('^(CalClustersAlg|second).*Energy')

  # accumulate the values
  raw_energy = 0
  corrected_energy = 0
  hashed_position = 0
  hashed_direction = 0
  file = open(filename)
  for line in file.readlines() :
    if regexp.search(line) :
      words = line.split()
      raw_energy += string.atof(words[3])
      corrected_energy += string.atof(words[5])
      hashed_position += string.atof(words[6]) + string.atof(words[7]) + string.atof(words[8])
      hashed_direction += string.atof(words[9]) + string.atof(words[9]) + string.atof(words[11])
      
  # result
  return [ raw_energy, corrected_energy, hashed_position, hashed_direction ]
   

#=================================================
# all things to be done for a given set of options
# prerequisite: the current dir is something
#   like <project>/<package>/<version>/cmt
#=================================================

def run_job(package,options,cmtbin_depend) :

  # change directory
  os.chdir(cmt_path(package))

  # file names
  exe_name = os.path.join('..',os.environ['CMTCONFIG'],'test_'+package+'.exe')
  opt_name = os.path.join(original_dir,options+'.txt')
  log_name = os.path.join(original_dir,options+'.log')
  ref_name = os.path.join(original_dir,options+'.ref')
  if cmtbin_depend == True :
    linux_expr = re.compile('^Linux')
    if linux_expr.search(os.environ['CMTBIN']) :
      ref_name = os.path.join(original_dir,options+'.linux.ref')
    windows_expr = re.compile('^VisualC')
    if windows_expr.search(os.environ['CMTBIN']) :
      ref_name = os.path.join(original_dir,options+'.win.ref')

  # command
  if os.name == 'posix':
    exe_command = '. setup.sh ; '+exe_name+' '+opt_name
  if os.name == 'nt':
    exe_command = 'call setup.bat & '+exe_name+' '+opt_name

  # prepare the log file
  log_file = file(log_name,'w')
  log_pipe = os.popen(exe_command)
  for line in log_pipe :
    log_file.write(line)
  log_file.close()
  if log_pipe.close() != None :
    print 'VALIDATION ERROR: '+package+' '+options+' EXECUTION FAILED'
    sys.exit(1)

  # compute energies
  log_key_values = key_values(log_name)
  ref_key_values = key_values(ref_name)

  # compare
  if log_key_values != ref_key_values :
    print 'VALIDATION ERROR: '+package+' '+options+' COMPARISON FAILED',log_key_values,'!=',ref_key_values
    #sys.exit(1)
  else :
    print 'validation: '+package+' '+options+' ok',log_key_values
    
  # back to original dir
  os.chdir(original_dir)


#=================================
#  jobs
#=================================

build_application_test('CalRecon')
build_application_test('Gleam')

run_job('CalRecon','jobOptions',False)
run_job('CalRecon','simpleOptions',False)
run_job('Gleam','gleamOptions',True)

