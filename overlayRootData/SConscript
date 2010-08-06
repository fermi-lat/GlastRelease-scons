# -*- python -*-
# $Header$
# Authors: Tracy Usher <usher@slac.stanford.edu>
# Version: overlayRootData-00-03-01
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

locIncs = listFiles(['overlayRootData/*.h'])
libEnv.Tool('addLinkDeps', package='overlayRootData', toBuild='rootlib')
overlayRootDataRootcint = libEnv.Rootcint('overlayRootData/overlayRootData_rootcint',
                                          ['overlayRootData/AcdOverlay.h',
                                           'overlayRootData/CalOverlay.h',
                                           'overlayRootData/EventOverlay.h',
                                           'overlayRootData/GlastAxis.h',
                                           'overlayRootData/PtOverlay.h',
                                           'overlayRootData/GemOverlay.h',
                                           'overlayRootData/TkrOverlay.h',
                                           'overlayRootData/LinkDef.h'],
                                          includes = [''],
                                          localIncludes = locIncs,
                                          packageName = 'overlayRootData')

libEnv['rootcint_node'] = overlayRootDataRootcint
                                       
overlayRootData = libEnv.RootDynamicLibrary('overlayRootData',
                                            listFiles(['src/*.cxx']) + ['overlayRootData/overlayRootData_rootcint.cxx'])

progEnv.Tool('overlayRootDataLib')
test_overlayRootData = progEnv.Program('test_overlayRootData',
                                       ['src/test/testOverlayClasses.cxx'])

progEnv.Tool('registerTargets', package = 'overlayRootData',
             rootcintSharedCxts = [[overlayRootData, libEnv]], 
             testAppCxts = [[test_overlayRootData, progEnv]],
             includes = listFiles(['overlayRootData/*.h']))




