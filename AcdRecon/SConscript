# -*- python -*-
# $Header$
# Authors: Heather Kelly <heather@slac.stanford.edu>, Eric Charles <echarles@slac.stanford.edu>
# Version: AcdRecon-04-03-00
import os
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('AcdReconLib', depsOnly = 1)
AcdReconLib = libEnv.SharedLibrary('AcdRecon', listFiles(['src/*.cxx', 'src/Dll/*.cxx']))

progEnv.Tool('AcdReconLib')
test_AcdRecon = progEnv.GaudiProgram('test_AcdRecon', listFiles(['src/test/*.cxx']), test = 1)

progEnv.Tool('registerObjects', package = 'facilities', libraries = [AcdReconLib], testApps = [test_AcdRecon],
             includes = listFiles(['AcdRecon/*.h']))
