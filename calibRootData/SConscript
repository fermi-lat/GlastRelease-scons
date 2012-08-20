# -*- python -*-
# $Header$
# Authors: Joanne Bogart <jrb@slac.stanford.edu>
# Version: calibRootData-01-08-00
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

locIncs = listFiles(['calibRootData/*'], recursive=True)
libEnv.Tool('addLinkDeps', package='calibRootData', toBuild='rootlib')
calibRootDataRootcint = libEnv.Rootcint('calibRootData/calibRootData_rootcint.cxx',
                                        ['calibRootData/DacCol.h',
                                         'calibRootData/Cal/CalDimension.h',
                                         'calibRootData/Cal/CalSerial.h',
                                         'calibRootData/Cal/CalGain.h',
                                         'calibRootData/Cal/CalGainCol.h',
                                         'calibRootData/Tkr/TkrTower.h',
                                         'calibRootData/Tkr/ChargeScale.h',
                                         'calibRootData/Tkr/Tot.h',
                                         'calibRootData/LinkDef.h'],
                                        includes = ['src/Cal', 'src/Acd',
                                                    'src/Tkr'],
                                        localIncludes=locIncs,
                                        packageName='calibRootData')
calibRootData = libEnv.RootDynamicLibrary('calibRootData',
                                          listFiles(['src/*.cxx',
                                                     'src/Cal/*.cxx',
                                                     'src/Tkr/*.cxx']) + ['calibRootData/calibRootData_rootcint.cxx'])

libEnv['rootcint_node'] = calibRootDataRootcint
progEnv.Tool('calibRootDataLib')

test_calibRootData = progEnv.Program('test_calibRootData',
                                     listFiles(['src/test/*.cxx']))
makeTot = progEnv.Program('makeTot',[ 'src/test/makeTot.cxx'])

progEnv.Tool('registerTargets', package = 'calibRootData',
             rootcintSharedCxts = [[calibRootData, libEnv]],
             testAppCxts = [[test_calibRootData, progEnv]],
             binaryCxts = [[makeTot, progEnv]],
             includes = listFiles(['calibRootData/*'], recursive = 1))




