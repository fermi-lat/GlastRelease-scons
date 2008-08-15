# -*- python -*-
# $Header$
# Authors: T. Burnett <tburnett@u.washington.edu>, David Chamont <chamont@poly.in2p3.fr>
# Version: digiRootData-11-10-00
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('digiRootDataLib', depsOnly = 1)


progEnv.Tool('digiRootDataLib')
ReadFilterStats = progEnv.Program('ReadFilterStats', ['apps/ReadFilterStats.cxx'])

progEnv.Tool('registerObjects', package = 'digiRootData', binaries = [ReadFilterStats])
