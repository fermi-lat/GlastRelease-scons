# -*- python -*-
# $Header$
# Authors:  Claudia Lavalley
# Version: gcrSelectRootData-02-02-05
import os
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

locIncs = listFiles(['gcrSelectRootData/*.h'])
libEnv.Tool('addLinkDeps', package='gcrSelectRootData', toBuild='rootlib')
gcrSelectRootDataRootcint=libEnv.Rootcint('gcrSelectRootData/gcrSelectRootData_rootcint',
                                          ['gcrSelectRootData/GcrSelectedXtal.h',
                                           'gcrSelectRootData/GcrSelectEvent.h',   
                                           'gcrSelectRootData/GcrSelect.h',   
                                           'gcrSelectRootData/GcrSelectVals.h',
                                           'gcrSelectRootData/LinkDef.h'],
                                          includes=[''],
                                          localIncludes = locIncs,
                                          packageName = 'gcrSelectRootData')
gcrSelectRootData = libEnv.RootDynamicLibrary('gcrSelectRootData',
                                              listFiles(['src/*.cxx']) + ['gcrSelectRootData/gcrSelectRootData_rootcint.cxx'])
progEnv.Tool('gcrSelectRootDataLib')
test_gcrSelectRootData  = progEnv.Program('test_gcrSelectRootData',
                                          ['src/test/testGcrClasses.cxx'])

progEnv.Tool('registerTargets', package = 'gcrSelectRootData',
             libraryCxts = [[gcrSelectRootData, libEnv]],
             testAppCxts = [[test_gcrSelectRootData, libEnv]],
             includes = listFiles(['gcrSelectRootData/*.h']))




