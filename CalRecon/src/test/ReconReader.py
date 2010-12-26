
import os
import sys
import ROOT
import numpy

from CalLayout import *

ROOT.gStyle.SetCanvasColor(ROOT.kWhite)

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
            for branch in self.MeritChain.GetListOfBranches():
                branchName = branch.GetName()
                branchTitle = branch.GetTitle()
                branchType = branchTitle.split('/')[-1].lower()
                if '[' not in branchName and branchType not in ['c']:
                    print 'Creating array for %s...' % branchTitle
                    a = numpy.array([0.0], dtype = branchType)
                    self.MeritChain.SetBranchAddress(branchName, a)
                    self.MeritArrayDict[branchName] = a
            self.MeritTreeFormula = ROOT.TTreeFormula('cut', cut,
                                                      self.MeritChain)
            print 'Done. %d entries found.' % self.MeritChain.GetEntries()
        
        # Relation file...
        print 'Opening the relation file...'
        if relFilePath is None:
            relFilePath = reconFilePath.replace('recon.root', 'rel.root')
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
        if i >= self.getEntries():
            sys.exit('No more events, bye!')
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

    def getCalCluster(self, i = 0):
        if i >= 0 and i < self.getNumClusters():
            return self.getCalClusterCol().At(i)

    def getCalUberCluster(self):
        return self.getCalCluster(self.getNumClusters() - 1)

    def getCalTotalNumXtals(self):
        return self.getCalXtalRecCol().GetEntries()

    def getCalClusterXtalList(self, i = 0):
        xtalList = []
        if self.RelTable is None:
            print 'No relation table available.'
            return xtalList
        cluster = self.getCalCluster(i)
        if cluster is None:
            print 'No CAL cluster at this time :-('
            return xtalList
        for relation in self.RelTable.getRelationTable():
            if isinstance(relation.getFirst(), ROOT.CalXtalRecData):
                xtal = relation.getFirst()
                if cluster == relation.getSecond():
                    xtalList.append(xtal)
        return xtalList


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
        if numClusters > 0:
            print 'Full info for the first cluster...'
            cluster = reader.getCalCluster(0)
            cluster.Print()
        answer = raw_input('Press q to quit, s to save or type a number...')
        try:
            eventNumber = int(answer)
        except:
            eventNumber += 1
