# -*- python -*-
# $Header$
# Authors: Heather Kelly <heather@lheapop.gsfc.nasa.gov>, Tracy Usher   <usher@slac.stanford.edu>
# Version: LdfConverter-04-03-02
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('LdfConverterLib', depsOnly = 1)
LdfConverter = libEnv.SharedLibrary('LdfConverter', listFiles(['src/*.cxx']) + listFiles(['src/Dll/*.cxx']))

progEnv.AppendUnique(CPPDEFINES = ['_FILE_OFFSET_BITS=64'])
progEnv.Tool('LdfConverterLib')

test_LdfConverter = progEnv.GaudiProgram('test_LdfConverter', ['src/test/*.cxx'], test = 1)

progEnv.Tool('registerObjects', package = 'LdfConverter', libraries = [LdfConverter], testApps = [test_LdfConverter], includes = listFiles(['LdfConverter/*.h']))
