# Attempt to CAL clustering using Minimal Spanning Trees.

import os
import sys
import ROOT
import math

from CalLayout   import *
from ReconReader import *



CAL_LAYOUT = CalLayout()
COLOR_LIST = [ROOT.kRed, ROOT.kBlue, ROOT.kGray + 1]
for i in range(100):
    COLOR_LIST.append(ROOT.kBlack)



class MSTNode:

    def __init__(self, xtal):
        self.Position = xtal.getPosition()
        self.Weight = xtal.getEnergy()
        self.MarkerXY = ROOT.TMarker(self.x(), self.y(), 20)
        self.MarkerXZ = ROOT.TMarker(self.x(), self.z(), 20)
        self.MarkerYZ = ROOT.TMarker(self.y(), self.z(), 20)
        self.setMarkerSize(0.75)

    def setMarkerColor(self, color):
        self.MarkerXY.SetMarkerColor(color)
        self.MarkerXZ.SetMarkerColor(color)
        self.MarkerYZ.SetMarkerColor(color)

    def setMarkerSize(self, size):
        self.MarkerXY.SetMarkerSize(size)
        self.MarkerXZ.SetMarkerSize(size)
        self.MarkerYZ.SetMarkerSize(size)

    def x(self):
        return self.Position.X()

    def y(self):
        return self.Position.Y()

    def z(self):
        return self.Position.Z()

    def w(self):
        return self.Weight

    def distToNode(self, other):
        return math.sqrt((self.x() - other.x())**2 +\
                         (self.y() - other.y())**2 +\
                         (self.z() - other.z())**2)

    def draw(self, view = 'xy'):
        if view == 'xy':
            self.MarkerXY.Draw()
        elif view == 'xz':
            self.MarkerXZ.Draw()
        elif view == 'yz':
            self.MarkerYZ.Draw()

    def __str__(self):
        return '(%.2f, %.2f, %.2f)' % (self.x(), self.y(), self.z())



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

    def getMSTNodes(self, nodeSet):
        minDist = 100000000000
        mstNode1 = None
        mstNode2 = None
        for node1 in self.NodeList:
            for node2 in nodeSet.NodeList:
                dist = node1.distToNode(node2)
                if dist < minDist:
                    minDist = dist
                    mstNode1 = node1
                    mstNode2 = node2
        return (mstNode1, mstNode2)



class MSTEdge:

    def __init__(self, node1, node2):
        self.Node1 = node1
        self.Node2 = node2
        self.Length = node1.distToNode(node2)
        self.LineXY = ROOT.TLine(node1.x(), node1.y(), node2.x(), node2.y())
        self.LineXZ = ROOT.TLine(node1.x(), node1.z(), node2.x(), node2.z())
        self.LineYZ = ROOT.TLine(node1.y(), node1.z(), node2.y(), node2.z())

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
        if self.Length > other.Length:
            return 1
        else:
            return -1

    def __str__(self):
        return '%s--%s, l = %.2f mm' % (self.Node1, self.Node2, self.Length)



class MST:

    def __init__(self):
        self.EdgeList = []
        self.NodeList = []
        self.WeightSum = 0.0
        self.__MeanEdgeLength = 0.0
        self.__RmsEdgeLength = 0.0

    def addEdge(self, edge):
        for node in [edge.Node1, edge.Node2]:
            self.addNode(node)
        self.EdgeList.append(edge)
        self.__MeanEdgeLength += edge.Length
        self.__RmsEdgeLength += (edge.Length)**2

    def getMeanEdgeLength(self):
        return self.__MeanEdgeLength/max(1, self.getNumEdges())

    def getRmsEdgeLength(self):
        numEdges = self.getNumEdges()
        if numEdges < 2:
            return 0
        else:
            var = self.__RmsEdgeLength/numEdges -\
                  (self.__MeanEdgeLength/numEdges)**2
            return math.sqrt(var)

    def addNode(self, node):
        if node not in self.NodeList:
            self.NodeList.append(node)
            self.WeightSum += node.w()

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

    def getMaxEdgeLength(self):
        longestEdge = self.getLongestEdge()
        if longestEdge is not None:
            return longestEdge.Length
        else:
            return 0.0

    def getNumEdges(self):
        return len(self.EdgeList)

    def getNumNodes(self):
        return len(self.NodeList)

    def setColor(self, color):
        for node in self.NodeList:
            node.setMarkerColor(color)
        for edge in self.EdgeList:
            edge.setLineColor(color)

    def setMarkerSize(self, size):
        for node in self.NodeList:
            node.setMarkerSize(size)

    def draw(self, view = 'xy'):
        for edge in self.EdgeList:
            edge.draw(view)
        for node in self.NodeList:
            node.draw(view)

    def __cmp__(self, other):
        if self.WeightSum > other.WeightSum:
            return 1
        else:
            return -1

    def __str__(self):
        return  '%d xtals, E = %.2f MeV, mean l = %.2f mm, rms l = %.2f mm' %\
               (self.getNumNodes(), self.WeightSum, self.getMeanEdgeLength(),
                self.getRmsEdgeLength())



class MSTClustering:

    def __init__(self, xtalCol, threshold):
        self.LengthThreshold = threshold
        self.TopCanvasUber = None
        self.TopCanvasClusters = None
        self.UberTree = MST()
        self.ClusterCol = []
        if xtalCol.GetEntries() == 0:
            print "No nodes found."
            return
        elif xtalCol.GetEntries() == 1:
            print 'Single node found, creating minimum spanning tree...'
            self.UberTree.addNode(MSTNode(xtalCol[0]))
        else:
            print 'Multiple nodes found, creating minimum spanning tree...'
            setA = MSTNodeSet([MSTNode(xtalCol[0])])
            setB = MSTNodeSet()
            for (i, xtal) in enumerate(xtalCol):
                if i != 0:
                    setB.addNode(MSTNode(xtal))
            while not setB.isEmpty():
                (node1, node2) = setA.getMSTNodes(setB)
                edge = MSTEdge(node1, node2)
                self.UberTree.addEdge(edge)
                if edge.Length > self.LengthThreshold:
                    edge.setLineStyle(7)
                setA.addNode(node2)
                setB.removeNode(node2)
        print 'Done, %d node(s) in the uber tree, total energy = %.2f MeV.' %\
              (self.getTotalNumNodes(), self.getTotalWeightSum())
        self.findClusters()

    def getTotalNumNodes(self):
        return self.UberTree.getNumNodes()

    def getTotalWeightSum(self):
        return self.UberTree.WeightSum

    def getNumClusters(self):
        return len(self.ClusterCol)

    def getCluster(self, index):
        return self.ClusterCol[index]

    def findClusters(self):
        print 'Clustering...'
        tree = MST()
        if self.getTotalNumNodes() == 1:
            tree.addNode(self.UberTree.NodeList[0])
        for edge in self.UberTree.getEdges():
            if edge.Length > self.LengthThreshold:
                tree.addNode(edge.Node1)
                self.ClusterCol.append(tree)
                tree = MST()
                tree.addNode(edge.Node2)
            else:
                tree.addEdge(edge)
        self.ClusterCol.append(tree)
        self.ClusterCol.sort(reverse = True)
        print 'Done, %d cluster(s) found.' % self.getNumClusters()

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
        self.TopCanvasClusters.cd()
        self.TopCanvasClusters.Clear()
        CAL_LAYOUT.draw('xy')
        for (i, cluster) in enumerate(self.ClusterCol):
            cluster.setColor(COLOR_LIST[i])
            cluster.draw('xy')
        self.TopCanvasClusters.Update()
            




if __name__ == '__main__':
    MAX_NUM_XTALS = 100
    LENGTH_THRESHOLD = 250
    #fp = '/data/mc/allGamma-GR-v18r4p2-OVRLY-L/skimPhilippe_OVRLY_recon.root'
    fp = '/data/mc/allGamma-GR-v18r4p2-FAKEOVRLY/skimPhilippe_FAKEOVRLY_recon.root'
    reader = ReconReader(fp)
    answer = ''
    evtNumber = 0
    while answer != 'q':
        reader.getEntry(evtNumber)
        xtalCol = reader.getCalXtalRecCol()
        numXtals = reader.getNumCalXtals()
        print '\nAnalyzing event %d, %d xtal(s) found.' % (evtNumber, numXtals)
        print 'Id = %d-%d, CalEnergyRaw = %.2f' %\
              (reader.getMeritVariable('EvtRun'),
               reader.getMeritVariable('EvtEventId'),
               reader.getMeritVariable('CalEnergyRaw'))
        if numXtals <= MAX_NUM_XTALS:
            clustering = MSTClustering(xtalCol, LENGTH_THRESHOLD)
            clustering.draw()
            for (i, c) in enumerate(clustering.ClusterCol):
                print '* Cluster %d: %s' % (i, c)
        else:
            print 'Too many xtals, skipping...'
        answer = raw_input('Press q to quit or type an event number...')
        try:
            evtNumber = int(answer)
        except:
            evtNumber += 1

