# -*- python -*-
# $Header$
# Authors: Joanne Bogart <jrb@slac.stanford.edu>
# Version: calibUtil-01-10-01
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('calibUtilLib', depsOnly = 1)
calibUtil = libEnv.SharedLibrary('calibUtil', listFiles(['src/*.cxx', 'src/*.c']))

progEnv.Tool('calibUtilLib')
test_strips = progEnv.Program('test_strips', ['src/test/test_strips.cxx', 'src/StripSrv.cxx', 'src/GenericSrv.cxx'])
test_meta = progEnv.Program('test_meta', 'src/test/test_meta.cxx')
calibCoverage = progEnv.Program('calibCoverage', ['src/dbIntegrity/calibCoverage.cxx', 'src/dbIntegrity/Coverage.cxx'])

progEnv.Tool('registerObjects', package = 'calibUtil', libraries = [calibUtil], testApps = [test_strips, test_meta], binaries = [calibCoverage], includes = listFiles(['calibUtil/*.h']))
