# -*- python -*-
# $Header$
# Authors: T.Burnett <tburnett@u.washington.edu>
# Version: Event-11-27-00
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('EventLib', depsOnly = 1)
EventLib = libEnv.SharedLibrary('Event', listFiles(['src/MonteCarlo/*.cxx', 'src/Recon/TkrRecon/*.cxx',
                                                    'src/Recon/CalRecon/*.cxx', 'src/Recon/AcdRecon/*.cxx',
                                                    'src/Digi/*.cxx', 'src/TopLevel/*.cpp', 'src/Utilities/*.cxx']))

progEnv.Tool('EventLib')
test_Event = progEnv.Program('test_Event', 'src/test/testmain.cxx')
test_Tables = progEnv.Program('test_Tables', 'src/test/test_RelTabs.cxx')
test_TkrRecon = progEnv.Program('test_TkrRecon', 'src/test/test_TkrRecon.cxx')

progEnv.Tool('registerObjects', package = 'Event', libraries = [EventLib], testApps = [test_Event, test_Tables, test_TkrRecon], includes = listFiles(['Event/*'], recursive = True))
