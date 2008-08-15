# -*- python -*-
# $Header$
# Authors: T. Burnett <tburnett@u.washington.edu>, R.Giannitrapani <riccardo@fisica.uniud.it>, Francesco Longo <Francesco.Longo@ts.infn.it>
# Version: G4Generator-05-22-00
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('G4GeneratorLib', depsOnly = 1)

G4Generator = libEnv.SharedLibrary('G4Generator', listFiles(['src/Utilities/*.cxx','src/Dll/*.cxx','src/*.cxx']))


progEnv.Tool('G4GeneratorLib')
test_G4Generator = progEnv.GaudiProgram('test_G4Generator', listFiles(['src/test/*.cxx']), test = 1)

progEnv.Tool('registerObjects', package = 'G4Generator', libraries = [G4Generator], testApps = [test_G4Generator], includes = listFiles(['G4Generator/*.h']))
