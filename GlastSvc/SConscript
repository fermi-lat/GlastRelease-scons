# -*- python -*-
# $Header$
# Authors: T. Burnett <tburnett@u.washington.edu>
# Version: GlastSvc-09-26-00
import os
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('GlastSvcLib', depsOnly = 1)
GlastSvcLib = libEnv.SharedLibrary('GlastSvc', listFiles(['src/Dll/*.cpp', 'src/GlastDetSvc/*.cxx',
                                                          'src/GlastRandomSvc/*.cxx', 'src/Reco/*.cxx',
                                                          'src/CurrentEventAlg/*.cxx', 'src/HistorySvc/*.cxx',
                                                          'src/EventSelector/*.cpp']))

progEnv.Tool('GlastSvcLib')
progEnv.Tool('identsLib')
test_GlastSvc = progEnv.GaudiProgram('test_GlastSvc', listFiles(['src/test/*.cxx']), test = 1)

progEnv.Tool('registerObjects', package = 'GlastSvc', libraries = [GlastSvcLib], testApps = [test_GlastSvc], includes = listFiles(['GlastSvc/*'], recursive = True))
