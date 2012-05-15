# -*- python -*-
# $Header$
# Authors:  Claudia Lavalley
# Version: gcrSelectRootData-02-02-08
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
libEnv['rootcint_node'] = gcrSelectRootDataRootcint

libsrcs = listFiles(['src/*.cxx'])
libsrcs.append('gcrSelectRootData/gcrSelectRootData_rootcint.cxx')

gcrSelectRootData = libEnv.RootDynamicLibrary('gcrSelectRootData', libsrcs)

progEnv.Tool('gcrSelectRootDataLib')
test_gcrSelectRootData  = progEnv.Program('test_gcrSelectRootData',
                                          ['src/test/testGcrClasses.cxx'])

progEnv.Tool('registerTargets', package = 'gcrSelectRootData',
             rootcintSharedCxts = [[gcrSelectRootData, libEnv]], 
             testAppCxts = [[test_gcrSelectRootData, progEnv]],
             includes = listFiles(['gcrSelectRootData/*.h']))




