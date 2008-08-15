# -*- python -*-
# $Header$ 
# Authors: David Chamont <chamont@llr.in2p3.fr>
# Version: rootUtil-01-03-01
import os
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('rootUtilLib', depsOnly = 1)

rootUtilRootcint = libEnv.Rootcint('rootUtil/rootUtil_rootcint', ['rootUtil/BranchGroup.h','rootUtil/BgDataHandleBase.h',
		'rootUtil/BgDataHandle.h', 'rootUtil/CompositeEventList.h','rootUtil/CelEventComponent.h','rootUtil/CelEventEntry.h',
  		'rootUtil/CelFileAndTreeSet.h','rootUtil/CelFileAndTreeOffset.h','rootUtil/CelEventIDs.h'])
 
rootUtil = libEnv.SharedLibrary('rootUtil', listFiles(['src/*.cxx']) + ['rootUtil/rootUtil_rootcint.cxx'])

progEnv.Tool('rootUtilLib')
testRootInterface = progEnv.Program('testRootInterface',['apps/testRootInterface.cxx'])
celRelocate = progEnv.Program('celRelocate',['apps/celRelocate.cxx'])
celConvert = progEnv.Program('celConvert',['apps/celConvert.cxx'])
celInspect = progEnv.Program('celInspect',['apps/celInspect.cxx'])
digiMeritRecon = progEnv.Program('digiMeritRecon',['apps/digiMeritRecon.cxx'])
mergeSkim = progEnv.Program('mergeSkim',['apps/mergeSkim.cxx'])
testGenerateData = progEnv.Program('testGenerateData',['apps/testGenerateData.cxx'])
testMakeCel = progEnv.Program('testMakeCel',['apps/testMakeCel.cxx'])
testReadCel = progEnv.Program('testReadCel',['apps/testReadCel.cxx'])


progEnv.Tool('registerObjects', package = 'rootUtil', libraries = [rootUtil], 
	testApps = [testMakeCel, testReadCel, testGenerateData, testRootInterface],
	binaries = [celRelocate,celConvert,celInspect,digiMeritRecon,mergeSkim],
	includes = listFiles(['rootUtil/*.h']))
