# Attempt to CAL clustering using Minimal Spanning Trees.

import os
import sys
import ROOT
import math
import time

from CalLayout          import *
from ReconReader        import *
from CalMomentsAnalysis import *
from ClusterPredictor   import *

#ROOT.gStyle.SetOptStat(111111)


CAL_LAYOUT = CalLayout()
COLOR_LIST = [ROOT.kRed, ROOT.kBlue, ROOT.kGray + 1]
for i in range(100):
    COLOR_LIST.append(ROOT.kBlack)



def nodeDist2(node1, node2):
    return (node1.X - node2.X)*(node1.X - node2.X) +\
           (node1.Y - node2.Y)*(node1.Y - node2.Y) +\
           (node1.Z - node2.Z)*(node1.Z - node2.Z)

def nodeDist(node1, node2):
    return math.sqrt(nodeDist2(node1, node2))

def getNextEdge(set1, set2):
    minWeight = 100000000000
    for n1 in set1.NodeList:
        for n2 in set2.NodeList:
            weight = nodeDist2(n1, n2)
            if weight < minWeight:
                minWeight = weight
                node1 = n1
                node2 = n2
    return MSTEdge(node1, node2, math.sqrt(minWeight))

def getWeigthThreshold(energy):
    if energy > 1000.0:
        return 200.0
    else:
        return 500.0 - 150.0*(math.log10(energy) - 1)





class MSTNode:

    def __init__(self, xtal):
        self.Position = xtal.getPosition()
        self.X = self.Position.X()
        self.Y = self.Position.Y()
        self.Z = self.Position.Z()
        self.Energy = xtal.getEnergy()
        self.PackedId = xtal.getPackedId()
        self.MarkerXY = ROOT.TMarker(self.X, self.Y, 20)
        self.MarkerXZ = ROOT.TMarker(self.X, self.Z, 20)
        self.MarkerYZ = ROOT.TMarker(self.Y, self.Z, 20)
        self.setMarkerSize(0.75)

    def setMarkerColor(self, color):
        self.MarkerXY.SetMarkerColor(color)
        self.MarkerXZ.SetMarkerColor(color)
        self.MarkerYZ.SetMarkerColor(color)

    def setMarkerSize(self, size):
        self.MarkerXY.SetMarkerSize(size)
        self.MarkerXZ.SetMarkerSize(size)
        self.MarkerYZ.SetMarkerSize(size)

    def draw(self, view = 'xy'):
        if view == 'xy':
            self.MarkerXY.Draw()
        elif view == 'xz':
            self.MarkerXZ.Draw()
        elif view == 'yz':
            self.MarkerYZ.Draw()

    # Basic interface to plug the moments analysis.

    def getPosition(self):
        return self.Position

    def getEnergy(self):
        return self.Energy

    def getPackedId(self):
        return self.PackedId

    def __str__(self):
        return '(%.2f, %.2f, %.2f)' % (self.X, self.Y, self.Z)



class MSTNodeSet:

    def __init__(self, nodeList = []):
        self.NodeList = []
        if len(nodeList):
            for node in nodeList:
                self.addNode(node)

    def getNumNodes(self):
        return len(self.NodeList)

    def getNode(self, i):
        return self.NodeList[i]

    def addNode(self, node):
        self.NodeList.append(node)

    def removeNode(self, node):
        self.NodeList.remove(node)

    def isEmpty(self):
        return self.getNumNodes() == 0



class MSTEdge:

    def __init__(self, node1, node2, weight):
        self.Node1 = node1
        self.Node2 = node2
        self.Weight = weight
        self.LineXY = ROOT.TLine(node1.X, node1.Y, node2.X, node2.Y)
        self.LineXZ = ROOT.TLine(node1.X, node1.Z, node2.X, node2.Z)
        self.LineYZ = ROOT.TLine(node1.Y, node1.Z, node2.Y, node2.Z)

    def setLineColor(self, color):
        self.LineXY.SetLineColor(color)
        self.LineXZ.SetLineColor(color)
        self.LineYZ.SetLineColor(color)

    def setLineWidth(self, width):
        self.LineXY.SetLineWidth(width)
        self.LineXZ.SetLineWidth(width)
        self.LineYZ.SetLineWidth(width)

    def setLineStyle(self, style):
        self.LineXY.SetLineStyle(style)
        self.LineXZ.SetLineStyle(style)
        self.LineYZ.SetLineStyle(style)

    def draw(self, view):
        if view == 'xy':
            self.LineXY.Draw()
        elif view == 'xz':
            self.LineXZ.Draw()
        elif view == 'yz':
            self.LineYZ.Draw()

    def __cmp__(self, other):
        if self.Weight > other.Weight:
            return 1
        else:
            return -1

    def __str__(self):
        return '%s--%s, l = %.2f mm' % (self.Node1, self.Node2, self.Weight)



class MST:

    def __init__(self, fakeReader = None):
        self.FakeReader = fakeReader
        self.EdgeList = []
        self.NodeList = []
        self.EnergySum = 0.0
        self.ProbDict = {}
        self.Color = ROOT.kBlack
        self.__MeanEdgeWeight = 0.0
        self.__RmsEdgeWeight = 0.0
        self.MomentsClusterInfo = None
        self.ScrewedUp = False
        self.OrigNumXtals = None

    def addEdge(self, edge):
        for node in [edge.Node1, edge.Node2]:
            self.addNode(node)
        self.EdgeList.append(edge)
        weight = edge.Weight
        self.__MeanEdgeWeight += weight
        self.__RmsEdgeWeight += weight*weight

    def getMeanEdgeWeight(self):
        return self.__MeanEdgeWeight/max(1, self.getNumEdges())

    def getRmsEdgeWeight(self):
        numEdges = self.getNumEdges()
        if numEdges < 2:
            return 0
        else:
            var = self.__RmsEdgeWeight/numEdges -\
                  (self.__MeanEdgeWeight/numEdges)**2
            return math.sqrt(var)

    def addNode(self, node):
        if node not in self.NodeList:
            self.NodeList.append(node)
            self.EnergySum += node.Energy

    def sortEdges(self):
        self.EdgeList.sort()

    def getEdge(self, index):
        return self.EdgeList[index]

    def getEdges(self):
        return self.EdgeList

    def getFirstEdge(self):
        return self.getEdge(0)

    def getLongestEdge(self):
        if self.getNumEdges() > 0:
            return max(self.EdgeList)
        else:
            return None

    def getMaxEdgeWeight(self):
        longestEdge = self.getLongestEdge()
        if longestEdge is not None:
            return longestEdge.Weight
        else:
            return 0.0

    def getNumEdges(self):
        return len(self.EdgeList)

    def getNumNodes(self):
        return len(self.NodeList)

    def setColor(self, color):
        self.Color = color
        for node in self.NodeList:
            node.setMarkerColor(color)
        for edge in self.EdgeList:
            edge.setLineColor(color)

    def setMarkerSize(self, size):
        for node in self.NodeList:
            node.setMarkerSize(size)

    def checkAgainstFake(self):
        if self.FakeReader is not None:
            self.OrigNumXtals = self.FakeReader.getMeritVariable('CalNumXtals')
            if abs(self.getNumNodes() - self.OrigNumXtals) > 3:
                self.ScrewedUp = True
                for edge in self.EdgeList:
                    edge.setLineStyle(3)

    def draw(self, view = 'xy'):
        for edge in self.EdgeList:
            edge.draw(view)
        for node in self.NodeList:
            node.draw(view)

    # Stuff related to the moments analysis

    def getTransRms(self):
        try:
            return math.sqrt(self.MomentsClusterInfo.RmsTrans/self.EnergySum)
        except AttributeError:
            return None

    def getLongRms(self):
        try:
            return math.sqrt(self.MomentsClusterInfo.RmsLong/self.EnergySum)
        except AttributeError:
            return None

    def getLongRmsAsym(self):
        try:
            return self.MomentsClusterInfo.RmsLongAsym
        except AttributeError:
            return None

    def getCentroid(self):
        try:
            return self.MomentsClusterInfo.Centroid
        except AttributeError:
            return None

    def getProb(self, topology):
        try:
            return self.ProbDict[topology]
        except KeyError:
            return -1

    def getProbsText(self):
        text = ''
        if self.ProbDict == {}:
            return text
        for topology in FILE_PATH_DICT.keys():
            text += 'P(%s) = %.3f\n' % (topology, self.getProb(topology))
        return text.strip('\n')

    def __cmp__(self, other):
        if self.FakeReader is not None and self.getCentroid() is not None and \
               other.getCentroid() is not None:
            x0 = self.FakeReader.getMeritVariable('CalXEcntr')
            y0 = self.FakeReader.getMeritVariable('CalYEcntr')
            z0 = self.FakeReader.getMeritVariable('CalZEcntr')
            x1 = self.getCentroid().X()
            y1 = self.getCentroid().Y()
            z1 = self.getCentroid().Z()
            x2 = other.getCentroid().X()
            y2 = other.getCentroid().Y()
            z2 = other.getCentroid().Z()
            dist1 = (x1 - x0)**2 + (y1 - y0)**2  + (z1 - z0)**2
            dist2 = (x2 - x0)**2 + (y2 - y0)**2  + (z2 - z0)**2
            if dist1 < dist2:
                return 1
            else:
                return -1
        if self.EnergySum > other.EnergySum:
            return 1
        else:
            return -1

    def __str__(self):
        text = '%d xtals, E = %.1f MeV, mean w = %.1f mm, rms w = %.1f mm' %\
               (self.getNumNodes(), self.EnergySum, self.getMeanEdgeWeight(),
                self.getRmsEdgeWeight())
        transRms = self.getTransRms()
        longRms = self.getLongRms()
        longRmsAsym = self.getLongRmsAsym()
        centroid = self.getCentroid()
        if transRms is not None:
            text += '\n  -> Trans. RMS = %.2f mm' % transRms
            text += ', long RMS = %.2f mm' % longRms
            text += ', long RMS asym = %.3f' % longRmsAsym
            text += '\n  -> Centroid at %s' % centroid
        return text



class MSTPredictor(ClusterPredictor):

    def __init__(self):
        ClusterPredictor.__init__(self, 'cluclass.root')

    def setCluster(self, cluster):
        self.Cluster = cluster

    def cutPassed(self):
        return True

    def getVar(self, varLabel):
        if varLabel == 'CalEnergyRaw':
            return self.Cluster.EnergySum
        elif varLabel == 'CalTransRms':
            return self.Cluster.getTransRms()
        elif varLabel == 'CalLRmsAsym':
            return self.Cluster.getLongRmsAsym()
        elif varLabel == 'MomentRatio':
            # 'log10(CalLongRms/CalTransRms)'
            try:
                return math.log10(self.Cluster.getLongRms()/\
                                  self.Cluster.getTransRms())
            except:
                return 0
        elif varLabel == 'NumXtals':
            # 'CalNumXtals/log10(CalEnergyRaw)'
            try:
                return self.Cluster.getNumNodes()/\
                       math.log10(self.Cluster.EnergySum)
            except:
                return 0

    def classifyEvent(self):
        try:
            ClusterPredictor.classifyEvent(self)
        except:
            print 'Could not classify event.'
            for topology in FILE_PATH_DICT.keys():
                self.ProbDict[topology] = 0.0


class MSTClustering:

    def __init__(self, xtalCol, threshold = None, doMomentAnalysis = True,
                 classifier = None, fakeReader = None):
        self.Classifier = classifier
        self.FakeReader = fakeReader
        startTime = time.time()
        self.WeightThreshold = threshold
        self.TopCanvasUber = None
        self.TopCanvasClusters = None
        self.EdgeWeightCanvas = None
        self.UberTree = MST(self.FakeReader)
        self.ClusterCol = []
        numXtals = xtalCol.GetEntries()
        if numXtals == 0:
            print "No nodes found."
            return
        elif numXtals == 1:
            print 'Single node found, that was easy.'
            self.UberTree.addNode(MSTNode(xtalCol[0]))
        else:
            print '%d nodes found, creating minimum spanning tree...' %\
                  numXtals
            setA = MSTNodeSet([MSTNode(xtalCol[0])])
            setB = MSTNodeSet()
            for (i, xtal) in enumerate(xtalCol):
                if i != 0:
                    setB.addNode(MSTNode(xtal))
            while not setB.isEmpty():
                edge = getNextEdge(setA, setB)
                self.UberTree.addEdge(edge)
                setA.addNode(edge.Node2)
                setB.removeNode(edge.Node2)
        elapsedTime = time.time() - startTime
        print 'Done in %.3f s, %d node(s) in the uber tree, E = %.1f MeV.' %\
              (elapsedTime, self.getTotalNumNodes(), self.getTotalEnergySum())
        self.findClusters(doMomentAnalysis)

    def getTotalNumNodes(self):
        return self.UberTree.getNumNodes()

    def getTotalEnergySum(self):
        return self.UberTree.EnergySum

    def getNumClusters(self):
        return len(self.ClusterCol)

    def getCluster(self, index):
        return self.ClusterCol[index]

    def getUberCluster(self):
        return self.UberTree

    def findClusters(self, doMomentAnalysis):
        startTime = time.time()
        self.WeightThreshold = self.WeightThreshold or \
                               getWeigthThreshold(self.getTotalEnergySum())
        print 'Doing clustering (weight threshold = %.2f mm)...' %\
              self.WeightThreshold
        tree = MST(self.FakeReader)
        if self.getTotalNumNodes() == 1:
            tree.addNode(self.UberTree.NodeList[0])
        for edge in self.UberTree.getEdges():
            if edge.Weight > self.WeightThreshold:
                edge.setLineStyle(7)
                tree.addNode(edge.Node1)
                self.ClusterCol.append(tree)
                tree = MST(self.FakeReader)
                tree.addNode(edge.Node2)
            else:
                tree.addEdge(edge)
        self.ClusterCol.append(tree)
        elapsedTime = time.time() - startTime
        print 'Done in %.3f s, %d cluster(s) found.' %\
              (elapsedTime, self.getNumClusters())
        if doMomentAnalysis:
            startTime = time.time()
            print 'Doing moment analysis...'
            for c in self.ClusterCol:
                c.MomentsClusterInfo = MomentsClusterInfo(c.NodeList)
            elapsedTime = time.time() - startTime
            print 'Done in %.3f s.' % elapsedTime
        if self.Classifier is not None:
            startTime = time.time()
            print 'Doing classification...'
            for c in self.ClusterCol:
                self.Classifier.setCluster(c)
                self.Classifier.classifyEvent()
                c.ProbDict = self.Classifier.ProbDict
            elapsedTime = time.time() - startTime
            print 'Done in %.3f s.' % elapsedTime
        self.ClusterCol.sort(reverse = True)    

    def draw(self):
        if self.TopCanvasUber is None:
            self.TopCanvasUber = getTopCanvas('Top (uber cluster)',
                                              topx = 10, topy = 10)
        if self.TopCanvasClusters is None:
            self.TopCanvasClusters = getTopCanvas('Top (all clusters)',
                                                  topx = 680, topy = 10)
        self.TopCanvasUber.cd()
        self.TopCanvasUber.Clear()
        CAL_LAYOUT.draw('xy')
        self.UberTree.draw('xy')
        self.TopCanvasUber.Update()
        self.TopCanvasUber.SaveAs('mst-uber.pdf')
        self.TopCanvasClusters.cd()
        self.TopCanvasClusters.Clear()
        CAL_LAYOUT.draw('xy')
        for (i, cluster) in enumerate(self.ClusterCol):
            cluster.setColor(COLOR_LIST[i])
            cluster.draw('xy')
        self.TextList = []
        for c in self.ClusterCol:
            try:
                x0 = c.getCentroid().X()
                y0 = c.getCentroid().Y()
                if x0 < 400:
                    x = x0 + 100
                else:
                    x = x0 - 400
                if y0 < 400:
                    y = y0 + 100
                else:
                    y = y0 - 200
                if c.ScrewedUp:
                    screwText = ' (were %d)' % c.OrigNumXtals
                else:
                    screwText = ''
                l = ROOT.TLatex(x, y + 50, '%d xtals%s, E = %d MeV' %\
                                (c.getNumNodes(), screwText, c.EnergySum))
                l.SetTextColor(c.Color)
                l.SetTextSize(15)
                l.Draw()
                self.TextList.append(l)
                for (i, topology) in enumerate(FILE_PATH_DICT.keys()):
                    l = ROOT.TLatex(x, y-50*i, 'P(%s) = %.2f' %\
                                    (topology, c.getProb(topology)))
                    l.SetTextColor(c.Color)
                    l.SetTextSize(15)
                    l.Draw()
                    self.TextList.append(l)
                #l = ROOT.TLatex(x, y, 'P(gam) = %.2f' % c.getProb('gam'))
                #l.SetTextColor(c.Color)
                #l.SetTextSize(15)
                #l.Draw()
                #self.TextList.append(l)
                #l = ROOT.TLatex(x, y - 50, 'P(had) = %.2f' % c.getProb('had'))
                #l.SetTextColor(c.Color)
                #l.SetTextSize(15)
                #l.Draw()
                #self.TextList.append(l)
                #l = ROOT.TLatex(x, y - 100, 'P(mip) = %.2f' % c.getProb('mip'))
                #l.SetTextColor(c.Color)
                #l.SetTextSize(15)
                #l.Draw()
                #self.TextList.append(l)
            except:
                pass
        self.TopCanvasClusters.Update()
        self.TopCanvasClusters.SaveAs('mst-clusters.pdf')




if __name__ == '__main__':
    from optparse import OptionParser
    usage = 'MSTClustering.py reconFile [opt] <reconFile> [fakeOvrlyFile]'
    parser = OptionParser()
    parser.add_option('-x', '--max-num-xtals', type = int, dest = 'x',
                      default = 10000,
                      help = 'maximum number of xtal for running the MST')
    parser.add_option('-w', '--max-edge-weight', type = float, dest = 'w',
                      default = None,
                      help = 'threshold length for the MST clustering (in mm)')
    parser.add_option('-M', '--moments-analysis', action = 'store_true',
                      dest = 'M', default = False,
                      help = 'run the moments analysis on the clusters')
    parser.add_option('-C', '--classification', action = 'store_true',
                      dest = 'C', default = False,
                      help = 'run the classification')
    (opts, args) = parser.parse_args()
    if len(args) == 0:
        print 'Usage: %s' % usage
        sys.exit('Please provide a recon input root file.')
    elif len(args) > 3:
        print 'Usage: %s' % usage
        sys.exit('Too many arguments.')
    inputFilePath = args[0]
    try:
        fakeOvrlyFile = args[1]
        fakeReader = ReconReader(fakeOvrlyFile)
    except IndexError:
        fakeOvrlyFile = None
        fakeReader = None
    reader = ReconReader(inputFilePath)
    if opts.C:
        classifier = MSTPredictor()
    else:
        classifier = None
    answer = ''
    evtNumber = 0
    while answer != 'q':
        reader.getEntry(evtNumber)
        if fakeReader is not None:
            fakeReader.getEntry(evtNumber)
        xtalCol = reader.getCalXtalRecCol()
        numXtals = reader.getNumCalXtals()
        print '\nAnalyzing event %d, %d xtal(s) found.' % (evtNumber, numXtals)
        runId = reader.getMeritVariable('EvtRun')
        evtId = reader.getMeritVariable('EvtEventId')
        print 'Id = %d-%d, CalEnergyRaw = %.2f' %\
              (runId, evtId, reader.getMeritVariable('CalEnergyRaw'))
        #longRms = reader.getMeritVariable('CalLongRms')
        #latRLn  = reader.getMeritVariable('CalLATRLn')
        #cntRLn  = reader.getMeritVariable('CalCntRLn')
        #logRLn  = math.log(latRLn - cntRLn)
        #print longRms, latRLn, cntRLn, longRms*logRLn
        if numXtals <= opts.x:
            clustering = MSTClustering(xtalCol, opts.w, opts.M, classifier,
                                       fakeReader)
            try:
                clustering.ClusterCol[0].checkAgainstFake()
            except IndexError:
                pass
            clustering.draw()
            for (i, c) in enumerate(clustering.ClusterCol):
                print '* Cluster %d: %s' % (i, c)
                print c.getProbsText()
        else:
            print 'Too many xtals (%d), skipping...' % opts.x
        answer = raw_input('Press q to quit, s to save or type a number...')
        try:
            evtNumber = int(answer)
        except:
            evtNumber += 1
        if answer == 's':
            filePath = 'mst_%d-%d_uber.pdf' % (runId, evtId)
            os.system('mv %s %s' % ('mst-uber.pdf', filePath))
            filePath = 'mst_%d-%d_clusters.pdf' % (runId, evtId)
            os.system('mv %s %s' % ('mst-clusters.pdf', filePath))
