# -*- python -*-
# $Header$
# Authors:  Claudia Lavalley
# Version: gcrSelectRootData-02-02-01
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

gcrSelectRootDataRootcint = libEnv.Rootcint('gcrSelectRootData/gcrSelectRootData_rootcint', ['gcrSelectRootData/GcrSelect.h',
                                                                                             'gcrSelectRootData/GcrSelectedXtal.h',
                                                                                             'gcrSelectRootData/GcrSelectEvent.h',
                                                                                             'gcrSelectRootData/GcrSelectVals.h'],
											     includes = [''])

libEnv.Tool('gcrSelectRootDataLib', depsOnly = 1)
gcrSelectRootData = libEnv.SharedLibrary('gcrSelectRootData', listFiles(['src/*.cxx']) + ['gcrSelectRootData/gcrSelectRootData_rootcint.cxx'])
progEnv.Tool('gcrSelectRootDataLib')
test_gcrSelectRootData  = progEnv.Program('test_gcrSelectRootData',['src/test/testGcrClasses.cxx'])

progEnv.Tool('registerObjects', package = 'gcrSelectRootData',
             libraries = [gcrSelectRootData],
             testApps = [test_gcrSelectRootData],
             includes = listFiles(['gcrSelectRootData/*.h']))




