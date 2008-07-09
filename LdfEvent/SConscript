# -*- python -*-
# $Header$
# Authors: Heather Kelly <heather@milkyway.gsfc.nasa.gov>
# Version: LdfEvent-04-05-00
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('LdfEventLib', depsOnly = 1)
LdfEvent = libEnv.SharedLibrary('LdfEvent', listFiles(['src/*.cxx']))

progEnv.Tool('LdfEventLib')
test_LdfEvent = progEnv.Program('test_LdfEvent', ['src/test/testMain.cxx'])

progEnv.Tool('registerObjects', package = 'LdfEvent', libraries = [LdfEvent], testApps = [test_LdfEvent], includes = listFiles(['LdfEvent/*.h']))
