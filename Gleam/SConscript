# -*- python -*-
# $Header$
# Authors: T. Burnett <tburnett@u.washington.edu>
# Version: Gleam-07-00-19

Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('GleamLib', depsOnly = 1)


progEnv.Tool('GleamLib')
Gleam = progEnv.Program('Gleam') #see corresponding requirements file.  not sure about that $(GlastMain) or how to use it.
test_Gleam = progEnv.GaudiProgram('test_Gleam', listFiles(['src/test/*.cxx']), test=1)
progEnv.Tool('registerObjects', package = 'Gleam', binaries = [Gleam], testApps = [test_Gleam], includes = listFiles(['Gleam/*.h']))
