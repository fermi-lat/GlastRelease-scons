
import array
import os

from MSTClustering import *


from optparse import OptionParser
usage = 'test_MST.py reconFile [opt] <reconFile>'
parser = OptionParser()
parser.add_option('-o', '--output-file', type = str, dest = 'o',
                  default = None,
                  help = 'path to the output file')
parser.add_option('-x', '--max-num-xtals', type = int, dest = 'x',
                  default = 10000,
                  help = 'maximum number of xtal for running the MST')
parser.add_option('-m', '--max-num-clusters', type = int, dest = 'm',
                  default = 5,
                  help = 'maximum number of clusters in the output tree')
parser.add_option('-w', '--max-edge-weight', type = float, dest = 'w',
                  default = None,
                  help = 'threshold length for the MST clustering (in mm)')
parser.add_option('-n', '--num-events', type = int, dest = 'n',
                  default = 10000000,
                  help = 'number of events to be processed')

(opts, args) = parser.parse_args()
if len(args) == 0:
    print 'Usage: %s' % usage
    sys.exit('Please provide a recon input root file.')
elif len(args) > 2:
    print 'Usage: %s' % usage
    sys.exit('Too many arguments.')
    
inputFilePath = args[0]
outputFilePath = opts.o or inputFilePath.replace('recon.root', 'MSTClu.root')
if os.path.exists(outputFilePath):
    sys.exit('Output file %s exists, remove it first.' % outputFilePath)
MAX_NUM_XTALS    = opts.x
WEIGHT_THRESHOLD = opts.w
MAX_NUM_CLUSTERS = opts.m

reader = ReconReader(inputFilePath)
numEvents = min(opts.n, reader.getEntries())
if numEvents < 0:
    numEvents = reader.getEntries()
outputFile = ROOT.TFile(outputFilePath, 'RECREATE')
outputTree = ROOT.TTree('MSTTuple', 'MSTTuple')
arrayDict = {}
BRANCH_DICT = {'EvtRun'             : ('i', 1),
               'EvtEventId'         : ('i', 1),
               'CalEnergyRaw'       : ('f', 1),
               'McEnergy'           : ('f', 1),
               'TkrNumTracks'       : ('f', 1),
               'CalCsIRLn'          : ('f', 1),
               'NumXtals'           : ('i', 1),
               'NumClusters'        : ('i', 1),
               'UberClusterNumXtals': ('f', 1),
               'UberClusterEnergy'  : ('f', 1),
               'UberClusterMeanW'   : ('f', 1),
               'UberClusterRmsW'    : ('f', 1),
               'UberClusterMaxW'    : ('f', 1),
               'ClusterNumXtals'    : ('i', MAX_NUM_CLUSTERS),
               'ClusterEnergy'      : ('f', MAX_NUM_CLUSTERS),
               'ClusterMeanW'       : ('f', MAX_NUM_CLUSTERS),
               'ClusterRmsW'        : ('f', MAX_NUM_CLUSTERS),
               'ClusterMaxW'        : ('f', MAX_NUM_CLUSTERS)
               }
for (branchName, (branchType, branchSize)) in BRANCH_DICT.items():
    a = array.array(branchType, [0]*branchSize)
    arrayDict[branchName] = a
    if branchSize == 1:
        branchTitle = '%s/%s' % (branchName, branchType.upper())
    else:
        branchTitle = '%s[%d]/%s' %\
                      (branchName, branchSize, branchType.upper())
    outputTree.Branch(branchName, a, branchTitle)
for i in xrange(numEvents):
    reader.getEntry(i)
    print '\nProcessing event %d/%d...' % (i, numEvents)
    xtalCol = reader.getCalXtalRecCol()
    numXtals = reader.getNumCalXtals()
    arrayDict['EvtRun'][0] = reader.getMeritVariable('EvtRun')
    arrayDict['EvtEventId'][0] = reader.getMeritVariable('EvtEventId')
    arrayDict['CalEnergyRaw'][0] = reader.getMeritVariable('CalEnergyRaw')
    arrayDict['McEnergy'][0] = reader.getMeritVariable('McEnergy')
    arrayDict['TkrNumTracks'][0] = reader.getMeritVariable('TkrNumTracks')
    arrayDict['CalCsIRLn'][0] = reader.getMeritVariable('CalCsIRLn')
    arrayDict['NumXtals'][0] = numXtals
    if numXtals <= MAX_NUM_XTALS:
        clustering = MSTClustering(xtalCol, WEIGHT_THRESHOLD)
        numClusters = clustering.getNumClusters()
        arrayDict['NumClusters'][0] = numClusters
        uberCluster = clustering.getUberCluster()
        arrayDict['UberClusterNumXtals'][0] = uberCluster.getNumNodes()
        arrayDict['UberClusterEnergy'][0] = uberCluster.EnergySum
        arrayDict['UberClusterMeanW'][0] = uberCluster.getMeanEdgeWeight()
        arrayDict['UberClusterRmsW'][0] = uberCluster.getRmsEdgeWeight()
        arrayDict['UberClusterMaxW'][0] = uberCluster.getMaxEdgeWeight()
        for cId in xrange(MAX_NUM_CLUSTERS):
            if cId < numClusters:
                c = clustering.getCluster(cId)
                arrayDict['ClusterNumXtals'][cId] = c.getNumNodes()
                arrayDict['ClusterEnergy'][cId] = c.EnergySum
                arrayDict['ClusterMeanW'][cId] = c.getMeanEdgeWeight()
                arrayDict['ClusterRmsW'][cId] = c.getRmsEdgeWeight()
                arrayDict['ClusterMaxW'][cId] = c.getMaxEdgeWeight()
            else:
                arrayDict['ClusterNumXtals'][cId] = 0
                arrayDict['ClusterEnergy'][cId] = 0.0
                arrayDict['ClusterMeanW'][cId] = 0.0
                arrayDict['ClusterRmsW'][cId] = 0.0
                arrayDict['ClusterMaxW'][cId] = 0.0
    else:
        arrayDict['NumClusters'][0] = 0
        arrayDict['UberClusterNumXtals'][0] = 0
        arrayDict['UberClusterEnergy'][0] = 0.0
        arrayDict['UberClusterMeanW'][0] = 0.0
        arrayDict['UberClusterRmsW'][0] = 0.0
        arrayDict['UberClusterMaxW'][0] = 0.0
        for cId in xrange(MAX_NUM_CLUSTERS):
            arrayDict['ClusterNumXtals'][cId] = 0
            arrayDict['ClusterEnergy'][cId] = 0.0
            arrayDict['ClusterMeanW'][cId] = 0.0
            arrayDict['ClusterRmsW'][cId] = 0.0
            arrayDict['ClusterMaxW'][cId] = 0.0
    outputTree.Fill()
outputFile.Write()
outputFile.Close()
