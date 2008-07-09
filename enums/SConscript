# -*- python -*-
# $Header$
# Authors: Heather Kelly <heather@milkyway.gsfc.nasa.gov>
# Version: enums-01-01-00
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()

test_enums = progEnv.Program('test_enums', 'src/test/test_enums.cxx')

progEnv.Tool('registerObjects', package = 'enums', testApps = [test_enums], includes = listFiles(['enums/*.h']))
