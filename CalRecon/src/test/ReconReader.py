
import os
import sys
import ROOT
import copy
import numpy

from CalLayout import *

from math import sqrt


ROOT.gStyle.SetCanvasColor(ROOT.kWhite)

LIBRARIES = ['libcommonRootData.so', 'libreconRootData.so']
MERIT_VARS = ['EvtEventId',
              'McEnergy', 'CalEnergyRaw', 'CTBBestEnergy', 'CalEnergyCorr',
              'McXDir', 'McYDir', 'McZDir',
              'Tkr1XDir', 'Tkr1YDir', 'Tkr1ZDir',
              'Tkr1X0', 'Tkr1Y0', 'Tkr1Z0',
              'CalCsIRLn', 'CalLATRLn',
              'CalTrackDoca', 'CalTrackAngle', 'CTBCORE',
              'CalTransRms', 'CalLongRms',
              'TkrNumTracks']

print 'Loading necessary libraries...'
LIB_DIRS = []
for var in ['RELEASE', 'PARENT']:
    if var in os.environ.keys():
        LIB_DIRS.append(os.path.join(os.environ[var], 'lib'))
if not len(LIB_DIRS):
    sys.exit('Could not locate libs, please define $RELEASE and/or $PARENT.')
for libName in LIBRARIES:
    for folder in LIB_DIRS:
        libPath = os.path.join(folder, libName)
        if os.path.exists(libPath):
            print 'Loading %s' % libPath
            ROOT.gSystem.Load(libPath)
            break
print 'Done.'




# Convenience class inheriting from ROOT.CalCluster.
# Beside wrapping some of the methods in order to save lines of code, it
# provides a __cmp__ method that allows to sort lists of objects.

class CalCluster(ROOT.CalCluster):

    def __init__(self, cluster):
        ROOT.CalCluster.__init__(self, cluster)

    def getEnergy(self):
        return self.getParams().getEnergy()

    def getCentroid(self):
        return self.getParams().getCentroid()

    def getAxis(self):
        return self.getParams().getAxis()

    def getTransRms(self):
        return sqrt(ROOT.CalCluster.getRmsTrans(self)/self.getEnergy())

    def __cmp__(self, other):
        return int(1000*(other.getEnergy() - self.getEnergy()))




class ReconReader:
    
    def __init__(self, reconFilePath, meritFilePath = None,
                 cut = ''):
        if not os.path.exists(reconFilePath):
            sys.exit('Could not find %s. Abort.' % reconFilePath)
        print 'Creating the recon chain...'
        self.ReconChain = ROOT.TChain('Recon')
        self.ReconChain.Add(reconFilePath)
        self.ReconEvent = ROOT.ReconEvent()
        self.ReconChain.SetBranchAddress('ReconEvent',\
                                         ROOT.AddressOf(self.ReconEvent))
        print 'Done. %d entries found.' % self.ReconChain.GetEntries()
        print 'Creating the merit chain...'
        self.MeritArrayDict = {}
        self.ExpectedSkewness = None
        if meritFilePath is None:
            meritFilePath = reconFilePath.replace('recon.root', 'merit.root')
            print 'No path specified, trying %s...' % meritFilePath
        if not os.path.exists(meritFilePath):
            print 'Could not find the merit, will ignore it.'
            self.MeritChain = None
        else:
            self.MeritChain = ROOT.TChain('MeritTuple')
            self.MeritChain.Add(meritFilePath)
            for branchName in MERIT_VARS:
                if branchName == 'EvtEventId':
                    a = numpy.array([0.0], dtype = 'l')
                else:
                    a = numpy.array([0.0], dtype = 'f')
                self.MeritChain.SetBranchAddress(branchName, a)
                self.MeritArrayDict[branchName] = a
            self.MeritTreeFormula = ROOT.TTreeFormula('cut', cut,
                                                      self.MeritChain)
            print 'Done. %d entries found.' % self.MeritChain.GetEntries()

    def getMeritVariable(self, branchName):
        try:
            return self.MeritArrayDict[branchName][0]
        except KeyError:
            return 'N/A'

    def getEventInfo(self):
        if self.MeritChain is None:
            info = 'Event info not available.'
        else:
            info = 'Event info:\n'
            for var in MERIT_VARS:
                info += '    %s = %s\n' % (var, self.getMeritVariable(var))
        info = info.strip('\n')
        return info

    def getEntry(self, i):
        print 'ReconReader retrieving event %d...' % i
        self.ReconChain.GetEvent(i)
        if self.MeritChain is not None:
            self.MeritChain.GetEntry(i)
            return self.MeritTreeFormula.EvalInstance()
        return 1

    def getCalRecon(self):
        return self.ReconEvent.getCalRecon()

    def getCalXtalRecCol(self):
        return self.getCalRecon().getCalXtalRecCol()

    def getNumClusters(self):
        return self.getCalClusterCol().GetEntries()

    def getCalClusterCol(self):
        return self.getCalRecon().getCalClusterCol()

    def getCalCluster(self, i):
        if i >= 0 and i < self.getNumClusters():
            return self.getCalClusterCol().At(i)

    def getCalUberCluster(self):
        return self.getCalCluster(self.getNumClusters() - 1)

    def getCalClusterParams(self, i):
        cluster = self.getCalCluster(i)
        if cluster is not None:
            return cluster.getParams()

    def getCalUberClusterParams(self):
        return self.getCalUberCluster().getParams()

    def getCalTotalNumXtals(self):
        return self.getCalXtalRecCol().GetEntries()

    def getCalClusterEnergy(self, i):
        params = self.getCalClusterParams(i)
        if params is not None:
            return params.getEnergy()

    def getCalClusterCentroid(self, i):
        params = self.getCalClusterParams(i)
        if params is not None:
            return params.getCentroid()

    def getCalClusterAxis(self, i):
        params = self.getCalClusterParams(i)
        if params is not None:
            return params.getAxis()

    def getCalClusterTransRms(self, i):
        cluster = self.getCalCluster(i)
        if cluster is not None:
            return sqrt(cluster.getRmsTrans()/self.getCalClusterEnergy(i))

    def getCalClusterNumTruncXtals(self, i):
        cluster = self.getCalCluster(i)
        if cluster is not None:
            return cluster.getNumTruncXtals()

    def getCalClusterNumSaturatedXtals(self, i):
        cluster = self.getCalCluster(i)
        if cluster is not None:
            return cluster.getNumSaturatedXtals()

    # This function is different from the others in that it returns a
    # (sorted!) python list of CalCluster objects rather than a TCollection
    # of ROOT.CalCluster objects.

    def getCalClusterList(self):
        list = [CalCluster(cluster) for cluster in self.getCalClusterCol()]
        list.sort()
        return list


    



if __name__ == '__main__':
    from optparse import OptionParser
    usage = 'ReconReader.py reconFile <meritFile>'
    parser = OptionParser()
    parser.add_option('-c', '--skim-cut', type = str, dest = 'c',
                      default = '1',
                      help = 'a cut to filter the events')
    (opts, args) = parser.parse_args()
    if len(args) == 0:
        print 'Usage: %s' % usage
        sys.exit('Please provide a recon input root file.')
    elif len(args) > 2:
        print 'Usage: %s' % usage
        sys.exit('Too many arguments.')
    reconFilePath = args[0]
    try:
        meritFilePath = args[1]
    except IndexError:
        meritFilePath = None
    fmt   = '%14s  '*5
    hline = '*'*79
    reader = ReconReader(reconFilePath, meritFilePath, opts.c)
    for event in xrange(10):
        reader.getEntry(event)
        #print reader.getEventInfo()
        numClusters = reader.getNumClusters()
        numXtals = reader.getCalTotalNumXtals()
        print '%d cluster(s) found, %d xtals in total.' %\
              (numClusters, numXtals)
        print hline
        print fmt % ('Cluster ID', 'Energy', 'Trans. rms', 'Trunc. xtals',
                     'Sat. xtals')
        print hline
        for (i, cluster) in enumerate(reader.getCalClusterList()):
            print fmt % (i, cluster.getEnergy(), cluster.getTransRms(),
                         cluster.getNumTruncXtals(),
                         cluster.getNumSaturatedXtals())
        print hline
