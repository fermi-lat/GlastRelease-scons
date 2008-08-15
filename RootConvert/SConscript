# -*- python -*-
# $Header$
# Authors: David Chamont <chamont@poly.in2p3.fr>, Heather Kelly <heather@slac.stanford.edu>
# Version: RootConvert-01-41-00
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('RootConvertLib', depsOnly = 1)

RootConvert = libEnv.SharedLibrary('RootConvert', listFiles(['src/Utilities/*.cxx','src/MonteCarlo/*.cxx',
							     'src/Digi/*.cxx','src/Recon/Acd*.cxx',
							     'src/Recon/Adf*.cxx','src/Recon/Cal*.cxx',
							     'src/Recon/Gcr*.cxx','src/GcrSelect/*.cxx']))

if baseEnv['PLATFORM'] == 'win32':
	progEnv.AppendUnique(CPPDEFINES = ['GLEAM', '__i386', 'EFC_FILTER', 'D_WIN32'])

progEnv.Tool('RootConvertLib')
test_RootConvert = progEnv.Program('test_RootConvert', ['src/test/test_RootConvert.cxx'])

progEnv.Tool('registerObjects', package = 'RootConvert', libraries = [RootConvert], testApps = [test_RootConvert], includes = listFiles(['RootConvert/*.h']))
