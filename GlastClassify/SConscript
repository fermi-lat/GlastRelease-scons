# -*- python -*-
# $Header$
# Authors: Tracy Usher <usher@slac.stanford.edu>
# Version: GlastClassify-06-05-01

Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('GlastClassifyLib', depsOnly = 1)

ClassificationTree = StaticLibrary('ClassificationTree', listFiles(['src/*.cxx'])  + listFiles(['src/ImActivityNodes/*.cxx']) + listFiles(['src/xmlBuilders/*.cxx']))


GlastClassify = libEnv.SharedLibrary('GlastClassify', listFiles(['src/Dll/*.cxx']))

progEnv.Tool('GlastClassifyLib')
apply = progEnv.Program('apply', listFiles(['src/apply/*.cxx']))
progEnv.Tool('registerObjects', package = 'GlastClassify', libraries = [ClassificationTree, GlastClassify], binaries = [apply], includes = listFiles(['GlastClassify/*.h']))
