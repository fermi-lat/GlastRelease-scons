# -*- python -*-
# $Header$
# Authors: H. Kelly <heather@slac.stanford.edu>, David Chamont <chamont@poly.in2p3.fr>
# Version: digiRootData-11-14-07
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

locIncs = listFiles(['digiRootData/*'])
libEnv.Tool('addLinkDeps', package='digiRootData', toBuild='rootlib')
digiRootDataRootcint = libEnv.Rootcint('digiRootData/digiRootData_rootcint',
                                       ['digiRootData/AcdDigi.h',
                                        'digiRootData/CalDiagnosticData.h',
                                        'digiRootData/CalDigi.h',
                                        'digiRootData/CalXtalReadout.h',
                                        'digiRootData/DigiEvent.h',
                                        'digiRootData/GlastAxis.h',
                                        'digiRootData/L1T.h',
                                        'digiRootData/GemTime.h',
                                        'digiRootData/TimeTone.h',
                                        'digiRootData/LsfTime.h',
                                        'digiRootData/RunInfo.h',
                                        'digiRootData/DatagramInfo.h',
                                        'digiRootData/GemScalers.h',
                                        'digiRootData/Configuration.h',
                                        'digiRootData/LsfKeys.h',
                                        'digiRootData/MetaEvent.h',
                                        'digiRootData/Ccsds.h',
                                        'digiRootData/TkrDiagnosticData.h',
                                        'digiRootData/TkrDigi.h',
                                        'digiRootData/Gem.h',
                                        'digiRootData/ErrorData.h',
                                        'digiRootData/Tem.h',
                                        'digiRootData/FilterStatus.h',
                                        'digiRootData/ObfFilterStatus.h',
                                        'digiRootData/AdfDigi.h',
                                        'digiRootData/EventSummaryData.h',
                                        'digiRootData/LinkDef.h'],
                                       includes = [''], localIncludes=locIncs,
                                       packageName='digiRootData')

libEnv['rootcint_node'] = digiRootDataRootcint
                                       
digiRootData = libEnv.RootDynamicLibrary('digiRootData',
                                         listFiles(['src/*.cxx']) + ['digiRootData/digiRootData_rootcint.cxx'])

progEnv.Tool('digiRootDataLib')
progEnv.Tool('facilitiesLib')
test_digiRootData = progEnv.Program('test_digiRootData',
                                    ['src/test/testDigiClasses.cxx'])
progEnv.Tool('addLibrary', library = baseEnv['ldfLibs'])
progEnv.Tool('lsfDataLib')
if baseEnv['PLATFORM'] == 'win32':
    XGetoptObj = '#/lib/' + baseEnv['VARIANT'] + '/XGetopt.obj'
    ReadFilterStats = progEnv.Program('ReadFilterStats',
                                      ['apps/ReadFilterStats.cxx',
                                       XGetoptObj])

else:
    ReadFilterStats = progEnv.Program('ReadFilterStats',
                                      ['apps/ReadFilterStats.cxx'])

progEnv.Tool('registerTargets', package = 'digiRootData',
             rootcintSharedCxts = [[digiRootData, libEnv]], 
             testAppCxts = [[test_digiRootData, progEnv]],
             binaryCxts=[[ReadFilterStats, progEnv]],
             includes = listFiles(['digiRootData/*.h']))




