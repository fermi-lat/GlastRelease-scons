# -*- python -*-
# $Header$
# Authors: Heather Kelly <heather@milkyway.gsfc.nasa.gov>,David Chamont <chamont@poly.in2p3.fr>
# Version: mcRootData-02-22-03
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

locIncs = listFiles(['mcRootData/*.h'])
libEnv.Tool('addLinkDeps', package='mcRootData', toBuild = 'rootlib')
mcRootDataRootcint = libEnv.Rootcint('mcRootData/mcRootData_rootcint',
                                     ['mcRootData/McIntegratingHit.h',
                                      'mcRootData/McParticle.h',
                                      'mcRootData/McPositionHit.h',
                                      'mcRootData/McEvent.h',
                                      'mcRootData/McTrajectory.h',
                                      'mcRootData/LinkDef.h'],
                                     includes = [''], localIncludes = locIncs,
                                     packageName = 'mcRootData')

mcRootData = libEnv.RootDynamicLibrary('mcRootData',
                                       listFiles(['src/*.cxx']) + ['mcRootData/mcRootData_rootcint.cxx'])
progEnv.Tool('mcRootDataLib')

test_mcRootData = progEnv.Program('test_mcRootData',
                                  ['src/test/testMcClasses.cxx'])

progEnv.Tool('registerTargets', package = 'mcRootData',
             rootcintSharedCxts = [[mcRootData, libEnv]],
             testAppCxts = [[test_mcRootData, progEnv]], 
             includes = listFiles(['mcRootData/*.h']))




