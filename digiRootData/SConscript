# -*- python -*-
# $Header$
# Authors: T. Burnett <tburnett@u.washington.edu>, David Chamont <chamont@poly.in2p3.fr>
# Version: digiRootData-11-12-00
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('digiRootDataLib', depsOnly = 1)
digiRootDataRootcint = libEnv.Rootcint('digiRootData/digiRootData_rootcint',['digiRootData/AcdDigi.h',
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
                                        'digiRootData/EventSummaryData.h', 'digiRootData/LinkDef.h'], includes = [''])

                                       
digiRootData = libEnv.SharedLibrary('digiRootData', listFiles(['src/*.cxx']) + ['digiRootData/digiRootData_rootcint.cxx'])

progEnv.Tool('digiRootDataLib')
test_digiRootData = progEnv.Program('test_digiRootData', ['src/test/testDigiClasses.cxx'])
progEnv.Tool('addLibrary', library = baseEnv['ldfLibs'])
progEnv.Tool('lsfDataLib')
ReadFilterStats = progEnv.Program('ReadFilterStats', ['apps/ReadFilterStats.cxx'])

progEnv.Tool('registerObjects', package = 'digiRootData', libraries = [digiRootData], 
	testApps = [test_digiRootData], binaries=[ReadFilterStats], includes = listFiles(['digiRootData/*.h']))




