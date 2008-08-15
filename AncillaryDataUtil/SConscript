# -*- python -*-
# $Header$
# Authors: N.Omodei <nicola.omodei@pi.infn.it>
# Version: AncillaryDataUtil-01-00-02
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('AncillaryDataUtilLib', depsOnly = 1)
AncillaryDataUtil = libEnv.StaticLibrary('AncillaryDataUtil', listFiles(['src/*.cxx']))

progEnv.Tool('AncillaryDataUtilLib')
test_AncillaryDataUtil = progEnv.GaudiProgram('test_AncillaryDataUtil',[[ 'src/test/testMain.cxx']], test = 1)

progEnv.Tool('registerObjects', package = 'AncillaryDataUtil', libraries = [AncillaryDataUtil], testApps = [test_AncillaryDataUtil], includes = listFiles(['AncillaryDataUtil/*.h']))

