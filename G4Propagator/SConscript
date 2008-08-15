# -*- python -*-
# $Header$
# Authors: Tracy Usher <usher@slac.stanford.edu>
# Version: G4Propagator-02-04-03

Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('G4PropagatorLib', depsOnly = 1)

G4Propagator = libEnv.SharedLibrary('G4Propagator', listFiles(['src/*.cxx']) + listFiles(['src/Dll/*.cxx']))

progEnv.Tool('G4PropagatorLib')
test_G4Propagator = progEnv.GaudiProgram('test_G4Propagator', listFiles(['src/test/*.cxx']), test=1)
progEnv.Tool('registerObjects', package = 'G4Propagator', libraries = [G4Propagator], testApps = [test_G4Propagator])
