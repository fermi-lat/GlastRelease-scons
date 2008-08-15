# -*- python -*-
# $Header$
# Authors: Tracy Usher <usher@SLAC.Stanford.edu>
# Version: OnboardFilterTds-00-09-00
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('OnboardFilterTdsLib', depsOnly = 1)

progEnv.AppendUnique(CPPDEFINES = ['GLEAM'])
progEnv.AppendUnique(CPPDEFINES = ['__i386'])
progEnv.AppendUnique(CPPDEFINES = ['EFC_FILTER'])

if baseEnv['PLATFORM'] == 'win32':
	progEnv.AppendUnique(CPPDEFINES = ['_WIN32'])
if baseEnv['PLATFORM'] != 'win32':
	progEnv.AppendUnique(CPPDEFINES = ['MSG_ST'])
	progEnv.AppendUnique(CPPDEFINES = ['__unix__'])
	progEnv.AppendUnique(CPPDEFINES = ['WUT_IMPLEMENTATION_PX'])	

OnboardFilterTds = libEnv.SharedLibrary('OnboardFilterTds', listFiles(['src/*.cxx']))

progEnv.Tool('OnboardFilterTdsLib')
test_OnboardFilterTds = progEnv.Program('test_OnboardFilterTds', listFiles(['src/test/*.cxx']))
progEnv.Tool('registerObjects', package = 'OnboardFilterTds', libraries = [OnboardFilterTds], testApps = [test_OnboardFilterTds], includes = listFiles(['OnboardFilterTds/*.h']))
