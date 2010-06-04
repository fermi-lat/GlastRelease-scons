
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
        self.XtalList = []

    def addXtal(self, xtal):
        self.XtalList.append(xtal)

    def getTotNumXtals(self):
        return len(self.XtalList)

    def getEnergy(self):
        return self.getParams().getEnergy()

    def getCentroid(self):
        return self.getParams().getCentroid()

    def getAxis(self):
        return self.getParams().getAxis()

    def getTransRms(self):
        return sqrt(ROOT.CalCluster.getRmsTrans(self)/self.getEnergy())

    def distToCentroid(self, cluster):
        diff = cluster.getCentroid()
        diff -= self.getCentroid()
        return diff.Mag()

    def distToAxis(self, cluster):
        diff = cluster.getCentroid()
        diff -= self.getCentroid()
        cross = cluster.getAxis().Cross(diff)
        return cross.Mag()        

    def __cmp__(self, other):
        if other.getEnergy() - self.getEnergy() > 0:
            return 1
        else:
            return -1




class ReconReader:
    
    def __init__(self, reconFilePath, meritFilePath = None,
                 relFilePath = None, cut = ''):
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
        print 'Opening the relation file...'
        if relFilePath is None:
            relFilePath = reconFilePath.replace('recon.root', 'relation.root')
            print 'No path specified, trying %s...' % relFilePath
        if not os.path.exists(relFilePath):
            print 'Could not find the relations, will ignore it.'
            self.RelationsChain = None
            self.RelTable = None
        else:
            self.RelationsChain = ROOT.TChain('Relations')
            self.RelationsChain.Add(relFilePath)
            self.RelTable = ROOT.RelTable()
            self.RelationsChain.SetBranchAddress('RelTable',\
                                                 ROOT.AddressOf(self.RelTable))
            print 'Done. %d entries found.' % self.RelationsChain.GetEntries()

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
        self.ReconChain.GetEvent(i)
        if self.RelationsChain is not None:
            self.RelationsChain.GetEvent(i)
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
    # If fillXtalList and the relations root file is available, the code
    # will dig into the relation table and fill the list fo xtals for each
    # cluster.

    def getCalClusterList(self, fillXtalList = True):
        if self.RelTable is None:
            fillXtalList = False
        clusterList = []
        if fillXtalList:
            clusterDict = {}
        for cluster in self.getCalClusterCol():
            calCluster = CalCluster(cluster)
            clusterList.append(calCluster)
            if fillXtalList:
                clusterDict[cluster] = calCluster
        clusterList.sort()
        if fillXtalList:
            for relation in self.RelTable.getRelationTable():
                if isinstance(relation.getKey(), ROOT.CalXtalRecData):
                    xtal = relation.getKey()
                    for cluster in relation.getValueCol():
                        calCluster = clusterDict[cluster]
                        calCluster.addXtal(xtal)
        return clusterList


        


    



if __name__ == '__main__':
    from optparse import OptionParser
    usage = 'ReconReader.py reconFile <meritFile> <relFile>'
    parser = OptionParser()
    parser.add_option('-c', '--skim-cut', type = str, dest = 'c',
                      default = '1',
                      help = 'a cut to filter the events')
    (opts, args) = parser.parse_args()
    if len(args) == 0:
        print 'Usage: %s' % usage
        sys.exit('Please provide a recon input root file.')
    elif len(args) > 3:
        print 'Usage: %s' % usage
        sys.exit('Too many arguments.')
    reconFilePath = args[0]
    try:
        meritFilePath = args[1]
    except IndexError:
        meritFilePath = None
    try:
        relFilePath = args[2]
    except IndexError:
        relFilePath = None
    fmt   = '%7s %12s %12s %13s %13s %13s'
    hline = '*'*79
    reader = ReconReader(reconFilePath, meritFilePath, relFilePath, opts.c)
    for event in xrange(10):
        print 'ReconReader retrieving event %d...' % event
        if reader.getEntry(event):
            #print reader.getEventInfo()
            numClusters = reader.getNumClusters()
            numXtals = reader.getCalTotalNumXtals()
            print '%d cluster(s) found, %d xtals in total.' %\
                  (numClusters, numXtals)
            print hline
            print fmt % ('Cluster', 'Energy', 'Trans. rms', 'Tot. xtals',
                         'Trunc. xtals', 'Sat. xtals')
            print hline
            clusterList = reader.getCalClusterList()
            uberCluster = clusterList[0]
            for (i, cluster) in enumerate(clusterList):
                print fmt % (i,
                             '%.3f' % cluster.getEnergy(),
                             '%.3f' % cluster.getTransRms(),
                             cluster.getTotNumXtals(),
                             cluster.getNumTruncXtals(),
                             cluster.getNumSaturatedXtals())
            print hline
    

