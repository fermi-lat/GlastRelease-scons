# -*- python -*-
# $Header$
# Authors: T.Burnett <tburnett@u.washington.edu>
# Version: geomrep-04-03-00
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('geomrepLib', depsOnly = 1)

geomrep = libEnv.StaticLibrary('geomrep', listFiles(['src/*.cxx']))

progEnv.Tool('geomrepLib')
progEnv.Tool('guiLib')
test_geomrep = progEnv.Program('testgeom', listFiles(['src/test/*.cxx']))

progEnv.Tool('registerObjects', package = 'geomrep', libraries = [geomrep], testApps = [test_geomrep], 
	includes = listFiles(['geomrep/*.h']))
