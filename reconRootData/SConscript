# -*- python -*-
# $Header$
# Authors: Heather Kelly <heather@milkyway.gsfc.nasa.gov>
# Version: reconRootData-09-25-01
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

locIncs = listFiles(['reconRootData/*.h'])
libEnv.Tool('addLinkDeps', package='reconRootData', toBuild='rootlib')
reconRootDataRootcint = libEnv.Rootcint('reconRootData/reconRootData_rootcint',
                                        ['reconRootData/AdfRecon.h',
                                         'reconRootData/AcdRecon.h',
                                         'reconRootData/AcdReconV2.h',
                                         'reconRootData/CalParams.h',
                                         'reconRootData/CalFitParams.h',
                                         'reconRootData/CalMSTreeParams.h',					 
                                         'reconRootData/CalClusterLayerData.h',
                                         'reconRootData/CalCluster.h',
                                         'reconRootData/CalCorToolResult.h',
                                         'reconRootData/CalEventEnergy.h',
                                         'reconRootData/CalRecon.h',
                                         'reconRootData/CalRangeRecData.h',
                                         'reconRootData/CalXtalRecData.h',
                                         'reconRootData/ReconEvent.h',
                                         'reconRootData/AcdEventTopology.h',
                                         'reconRootData/AcdTkrIntersection.h',
                                         'reconRootData/AcdTkrPoca.h',
                                         'reconRootData/AcdHit.h',
                                         'reconRootData/AcdPocaData.h',
                                         'reconRootData/AcdTkrLocalCoords.h',
                                         'reconRootData/AcdTkrPoint.h',
                                         'reconRootData/AcdTkrHitPoca.h',
                                         'reconRootData/AcdTkrGapPoca.h',
                                         'reconRootData/AcdSplashVars.h',
                                         'reconRootData/AcdPocaDataV2.h', 
                                         'reconRootData/AcdTkrLocalCoordsV2.h',
                                         'reconRootData/AcdTkrPointV2.h',
                                         'reconRootData/AcdTkrHitPocaV2.h',
                                         'reconRootData/AcdTkrGapPocaV2.h',
                                         'reconRootData/TkrCluster.h',
                                         'reconRootData/TkrTrack.h',
                                         'reconRootData/TkrTrackHit.h',
                                         'reconRootData/TkrTrackParams.h',
                                         'reconRootData/TkrRecon.h',
                                         'reconRootData/TkrVertex.h',
                                         'reconRootData/TkrTruncationData.h',
                                         'reconRootData/TkrDiagnostics.h',
                                         'reconRootData/TaggerCluster.h',
                                         'reconRootData/AdfRecon.h',
                                         'reconRootData/GcrXtal.h',
                                         'reconRootData/GcrTrack.h',
					 'reconRootData/LinkDef.h'], 
                                        includes = [''],
                                        localIncludes = locIncs,
                                        packageName = "reconRootData")

libEnv['rootcint_node'] = reconRootDataRootcint
reconRootData = libEnv.SharedLibrary('reconRootData',
                                     listFiles(['src/*.cxx']) + ['reconRootData/reconRootData_rootcint.cxx'])
progEnv.Tool('reconRootDataLib')

test_reconRootData = progEnv.Program('test_reconRootData',
                                     ['src/test/testReconClasses.cxx'])
progEnv.Tool('registerTargets', package = 'reconRootData',
             rootcintSharedCxts = [[reconRootData, libEnv]], 
             testAppCxts = [[test_reconRootData, progEnv]],
             includes = listFiles(['reconRootData/*.h']))




