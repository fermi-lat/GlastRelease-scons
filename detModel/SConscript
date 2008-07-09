# -*- python -*-
# $Header$
# Authors: Riccardo Giannitrapani <riccardo@fisica.uniud.it>, Joanne Bogart <jrb@slac.stanford.edu>
# Version: detModel-02-21-00
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('detModelLib', depsOnly = 1)
detModel = libEnv.StaticLibrary('detModel', listFiles(['src/Management/*.cxx', 'src/Sections/*.cxx', 'src/Constants/*.cxx',
                                                       'src/Utilities/*.cxx', 'src/Materials/*cxx', 'src/Gdd.cxx']))

progEnv.Tool('detModelLib')
progEnv.Tool('xmlUtilLib')
progEnv.Tool('identsLib')
test_detModel = progEnv.Program('test_detModel', 'src/main.cxx')

progEnv.Tool('registerObjects', package = 'detModel', libraries = [detModel], testApps = [test_detModel], includes = listFiles(['detModel/*'], recursive = 1))
