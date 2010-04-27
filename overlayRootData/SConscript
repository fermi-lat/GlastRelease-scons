# -*- python -*-
# $Header$
# Authors: Tracy Usher <usher@slac.stanford.edu>
# Version: overlayRootData-00-03-00
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('overlayRootDataLib', depsOnly = 1)
overlayRootDataRootcint = libEnv.Rootcint('overlayRootData/overlayRootData_rootcint',
                                          ['overlayRootData/AcdOverlay.h',
                                           'overlayRootData/CalOverlay.h',
                                           'overlayRootData/EventOverlay.h',
                                           'overlayRootData/GlastAxis.h',
                                           'overlayRootData/PtOverlay.h',
                                           'overlayRootData/GemOverlay.h',
                                           'overlayRootData/TkrOverlay.h',
                                           'overlayRootData/LinkDef.h'],
                                          includes = [''])

libEnv['rootcint_node'] = overlayRootDataRootcint
                                       
overlayRootData = libEnv.SharedLibrary('overlayRootData',
                                       listFiles(['src/*.cxx']) + ['overlayRootData/overlayRootData_rootcint.cxx'])

progEnv.Tool('overlayRootDataLib')
test_overlayRootData = progEnv.Program('test_overlayRootData',
                                       ['src/test/testOverlayClasses.cxx'])

progEnv.Tool('registerTargets', package = 'overlayRootData',
             rootcintSharedCxts = [[overlayRootData, libEnv]], 
             testAppCxts = [[test_overlayRootData, progEnv]],
             includes = listFiles(['overlayRootData/*.h']))




