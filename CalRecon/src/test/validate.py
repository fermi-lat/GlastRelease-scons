#!/usr/bin/env python
import sys, string, os, re ;

# Depends on Gleam package !

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

  # check if it is the original package
 # if original_dir.rfind(os.path.join('',package,'')) != -1  :
 #   return os.path.join('..','..')
  
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
  print 'VALIDATION ERROR: package',package,'NOT FOUND'
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
  regexp = re.compile('^(test_CalRecon).*Energy')

  # accumulate the values
  nb_clusters = 0
  nb_profiles = 0
  nb_ll = 0
  raw_energy = 0
  corrected_energy = 0
  ll_energy = 0
  fit_energy = 0
  hashed_posdir = 0
  file = open(filename)
  for line in file.readlines() :
    if regexp.search(line) :
      words = line.split()
      if words[2] == 'Profile' :
        if words[5] != '0' :
          nb_profiles += 1
          fit_energy += string.atof(words[5])
      elif words[2] == 'Last' :
        energy = string.atof(words[6])
        if energy > 0 :
          nb_ll += 1
          ll_energy += energy
      elif words[3] != '0' :
        nb_clusters += 1
        raw_energy += string.atof(words[3])
        corrected_energy += string.atof(words[5])
        hashed_posdir += string.atof(words[6]) + string.atof(words[7]) + string.atof(words[8])
        hashed_posdir += string.atof(words[9]) + string.atof(words[9]) + string.atof(words[11])
      
  # result
  if nb_clusters>0 :
    raw_energy = raw_energy/nb_clusters
    corrected_energy = corrected_energy/nb_clusters
    hashed_posdir = hashed_posdir/nb_clusters
  if nb_profiles>0 :
    fit_energy = fit_energy/nb_profiles
  if nb_ll>0 :
    ll_energy = ll_energy/nb_ll
  ft_count = '%d'
  ft_mean = '%g'
  all_values = [ (ft_count % nb_clusters)+'*'+ (ft_mean % raw_energy) ]
  all_values = all_values + [ (ft_count % nb_ll)+'*'+ (ft_mean % ll_energy) ]
  all_values = all_values + [ (ft_count % nb_profiles)+'*'+ (ft_mean % fit_energy) ]
  all_values = all_values + [ (ft_count % nb_clusters)+'*'+ (ft_mean % corrected_energy) ]
  all_values = all_values + [ ft_mean % hashed_posdir ]
  return all_values

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
    print 'VALIDATION ERROR: '+binary_package+' '+options+' EXECUTION FAILED'
    sys.exit(1)

  # compute energies
  log_key_values = key_values(log_name)
  ref_key_values = key_values(ref_name)

  # compare
  if log_key_values != ref_key_values :
    print 'VALIDATION ERROR: '+binary_package+' '+options+' COMPARISON FAILED\n  ',log_key_values,'!=\n  ',ref_key_values
    #sys.exit(1)
  else :
    print 'validation: '+binary_package+' '+options+' ok',log_key_values
    
  # back to original dir
  os.chdir(original_dir)


#=================================
#  jobs
#=================================

# build the application
build_application_test('CalRecon')

# pure CaloRecon jobs 
run_job('CalRecon','CalRecon','jobOptions',False)
run_job('CalRecon','CalRecon','simpleOptions',False)

