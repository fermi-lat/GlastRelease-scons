# -*- python -*-
# $Header$
# Authors: M. Kocian <kocian@slac.stanford.edu>
# Version: configData-01-10-04
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('configDataLib', depsOnly = 1)

configData = libEnv.SharedLibrary('configData', listFiles(['configData/*.cxx']) + listFiles(['src/gem/*.cxx']) + listFiles(['src/fsw/*.cxx']) + \
listFiles(['src/base/*.cxx']) + listFiles(['src/fsw/*.db']))


progEnv.Tool('configDataLib')

if baseEnv['PLATFORM'] != 'win32':
	progEnv.AppendUnique(CCFLAGS = ['-Wno-unused-parameter'])



ConfigXml2Root = progEnv.Program('ConfigXml2Root',['apps/ConfigXml2Root.cxx'])
RootConfigCompare = progEnv.Program('RootConfigCompare',['apps/RootConfigCompare.cxx'])
findConfigs = progEnv.Program('findConfigs',['apps/findConfigs.cxx'])
dumpGemConfiguration = progEnv.Program('dumpGemConfiguration',['apps/dumpGemConfiguration.cxx'])
testdb = progEnv.Program('testdb',['apps/testdb.cxx'])
test_configData = progEnv.Program('test_configData',['src/test/test_configData.cxx'])

progEnv.Tool('registerObjects', package = 'configData', libraries = [configData], testApps = [test_configData], binaries = [ConfigXml2Root,RootConfigCompare,findConfigs \
,dumpGemConfiguration,testdb])
