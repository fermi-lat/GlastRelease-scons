#!/usr/bin/env python
import sys, string, os, re ;

# This script depends on Gleam package
# It must be run whenever the generation
# process or the digis format change.
# The result is gamma_5_gev_normal.root,
# which is used by validate.py.

#=================================================
# globals
#=================================================

original_dir = os.getcwd()
release_expr = re.compile('^v[0-9]+r[0-9]+')


#=================================================
# the function below establishes the path
# to the root directory of a given client package
#=================================================

def root_dir(package) :

  # ask cmt
  packages_pipe = os.popen('cd '+original_dir+' ; cmt show packages')
  for line in packages_pipe :
    tokens = line.split()
    if tokens[0] == package :
      if tokens[1] == 'v1' :
        packages_pipe.close()
        return os.path.join(tokens[2],tokens[0])
      else :
        packages_pipe.close()
        return os.path.join(tokens[2],tokens[0],tokens[1])
  packages_pipe.close()

  # not found
  print 'PREPARATION ERROR: package',package,'NOT FOUND'
  sys.exit(1)
      

#=================================================
# build a package test application
#=================================================

def build_application_test(package) :

  os.chdir(os.path.join(root_dir(package),'cmt'))
  
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
#  elif os.name == 'nt':
#    build_command += ' & cmt bro -local nmake /f nmake'
#    build_command += ' & nmake /f nmake test'
#    if os.system(build_command) != 0 :
#      print 'VALIDATION ERROR: test_'+package+'.exe BUILD FAILED'
#      sys.exit(1)

  os.chdir(original_dir)
   

#=================================================
# all things to be done for a given set of options
# prerequisite: the current dir is something
#   like <project>/<package>/<version>/cmt
#=================================================

def run_job(setup_package,binary_package,options,cmtbin_depend) :

  # change directory
  os.chdir(os.path.join(root_dir(setup_package),'cmt'))

  # file names
  exe_name = os.path.join(root_dir(binary_package),os.environ['CMTCONFIG'],'test_'+binary_package+'.exe')
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
    print 'PREPARATION ERROR: '+binary_package+' '+options+' EXECUTION FAILED'
    sys.exit(1)

  # back to original dir
  os.chdir(original_dir)


#=================================
#  job
#=================================

build_application_test('CalRecon')
run_job('Gleam','CalRecon','gamma_5_gev_normal',True)

