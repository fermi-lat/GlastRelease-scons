
import os
import sys
import ROOT
import copy
import numpy

from CalLayout import *

from math import sqrt, acos, log10


ROOT.gStyle.SetCanvasColor(ROOT.kWhite)

MERIT_VARS = ['EvtEventId', 'EvtRun',
              'McEnergy', 'CalEnergyRaw', 'CTBBestEnergy', 'CalEnergyCorr',
              'McXDir', 'McYDir', 'McZDir',
              'Tkr1XDir', 'Tkr1YDir', 'Tkr1ZDir',
              'CalXDir', 'CalYDir', 'CalZDir',
              'CalXEcntr', 'CalYEcntr', 'CalZEcntr',
              'CalUberXDir', 'CalUberYDir', 'CalUberZDir',
              'Tkr1X0', 'Tkr1Y0', 'Tkr1Z0',
              'CalCsIRLn', 'CalLATRLn', 'CalCntRLn',
              'CalTrackDoca', 'CalTrackAngle', 'CTBCORE',
              'CalTransRms', 'CalLongRms', 'CalLRmsAsym', 'CalNumXtals',
              'TkrNumTracks', 'TkrEnergy']

#
# Load the recon libraries.
#
LIBRARIES = ['libcommonRootData.so', 'libreconRootData.so']
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

#
# Python interface to read a recon file and have access to the objects
# in there. It has the capability of parsing the content of the
# corresponding merit and relation files, if passed as arguments.
#

class ReconReader:
    
    def __init__(self, reconFilePath, meritFilePath = None,
                 relFilePath = None, cut = '1'):
        # Recon file...
        if not os.path.exists(reconFilePath):
            sys.exit('Could not find %s. Abort.' % reconFilePath)
        print 'Creating the recon chain...'
        self.ReconChain = ROOT.TChain('Recon')
        self.ReconChain.Add(reconFilePath)
        self.ReconEvent = ROOT.ReconEvent()
        reconAddress = ROOT.AddressOf(self.ReconEvent)
        self.ReconChain.SetBranchAddress('ReconEvent', reconAddress)
        print 'Done. %d entries found.' % self.ReconChain.GetEntries()

        # Merit file...
        print 'Creating the merit chain...'
        self.MeritArrayDict = {}
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
                if branchName in ['EvtEventId', 'EvtRun']:
                    a = numpy.array([0.0], dtype = 'l')
                else:
                    a = numpy.array([0.0], dtype = 'f')
                self.MeritChain.SetBranchAddress(branchName, a)
                self.MeritArrayDict[branchName] = a
            self.MeritTreeFormula = ROOT.TTreeFormula('cut', cut,
                                                      self.MeritChain)
            print 'Done. %d entries found.' % self.MeritChain.GetEntries()
        
        # Relation file...
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
            relAddress = ROOT.AddressOf(self.RelTable)
            self.RelationsChain.SetBranchAddress('RelTable', relAddress)
            print 'Done. %d entries found.' % self.RelationsChain.GetEntries()

    def getEntries(self):
        return self.ReconChain.GetEntries()

    def getEntry(self, i):
        self.ReconChain.GetEvent(i)
        if self.RelationsChain is not None:
            self.RelationsChain.GetEvent(i)
        if self.MeritChain is not None:
            self.MeritChain.GetEntry(i)
            return self.MeritTreeFormula.EvalInstance()
        return 1

    def getMeritVariable(self, branchName):
        try:
            return self.MeritArrayDict[branchName][0]
        except KeyError:
            return None

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

    def getCalTotalNumXtals(self):
        return self.getCalXtalRecCol().GetEntries()

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
    reader = ReconReader(reconFilePath, meritFilePath, relFilePath, opts.c)
    answer = ''
    eventNumber = 0
    while answer != 'q':
        print 'ReconReader retrieving event %d...' % eventNumber
        reader.getEntry(eventNumber)
        numClusters = reader.getNumClusters()
        numXtals = reader.getCalTotalNumXtals()
        print '%d cluster(s) found, %d xtals in total.' %\
            (numClusters, numXtals)
        print 'Full info for the first cluster...'
        reader.getCalCluster(0).Print()
        answer = raw_input('Press q to quit, s to save or type a number...')
        try:
            eventNumber = int(answer)
        except:
            eventNumber += 1
    

"""
----------------------------------------------------
----------- Output from the fitting tool -----------
----------------------------------------------------
Energy = -1 +- -1 MeV
Centroid = (-330.76, 618.478, -180.831) mm
Centroid covariance matrix:
| 1  0  0 |
| 0  1  0 |
| 0  0  1 |
Axis = (0.362045, 0.686613, 0.630465)
Axis covariance matrix:
| 1  0  0 |
| 0  1  0 |
| 0  0  1 |Number of layers for the fit: 4
Fit chisquare: 1.12321e-06
"""
