# -*- python -*-
# $Header$
# Authors:  Claudia Lavalley
# Version: gcrSelectRootData-02-02-04
import os
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()
libEnv.Tool('gcrSelectRootDataLib', depsOnly = 1)
gcrSelectRootDataRootcint=libEnv.Rootcint('gcrSelectRootData/gcrSelectRootData_rootcint',
                                          ['gcrSelectRootData/GcrSelectedXtal.h',
                                           'gcrSelectRootData/GcrSelectEvent.h',   
                                           'gcrSelectRootData/GcrSelect.h',   
                                           'gcrSelectRootData/GcrSelectVals.h',
                                           'gcrSelectRootData/LinkDef.h'],
                                          includes=[''])
gcrSelectRootData = libEnv.SharedLibrary('gcrSelectRootData',
                                         listFiles(['src/*.cxx']) + ['gcrSelectRootData/gcrSelectRootData_rootcint.cxx'])
progEnv.Tool('gcrSelectRootDataLib')
test_gcrSelectRootData  = progEnv.Program('test_gcrSelectRootData',
                                          ['src/test/testGcrClasses.cxx'])

progEnv.Tool('registerTargets', package = 'gcrSelectRootData',
             libraryCxts = [[gcrSelectRootData, libEnv]],
             testAppCxts = [[test_gcrSelectRootData, libEnv]],
             includes = listFiles(['gcrSelectRootData/*.h']))




