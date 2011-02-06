
import math

from ReconReader import *
from CalLayout   import *

CSI_RAD_LEN = 18.6
CSI_MOL_RAD = 35.7

CAL_TOWER_GAP = TOWER_PITCH - 0.5*(CSI_WIDTH*12 + CSI_LENGTH)

COLOR_WHEEL = [ROOT.kRed, ROOT.kBlue, ROOT.kMagenta, ROOT.kGreen]


class Point(ROOT.TVector3):

    def __init__(self, x = 0., y = 0., z = 0.):
        ROOT.TVector3.__init__(self, x, y, z)
    
    def x(self):
        return self.X()

    def y(self):
        return self.Y()

    def z(self):
        return self.Z()

    def mag2(self):
        return self.Mag2()

    def __sub__(self, other):
        return Point(self.x() - other.x(),
                     self.y() - other.y(),
                     self.z() - other.z())

    def __mul__(self, other):
        return Point(self.x()*other, self.y()*other, self.z()*other)

    def __neg__(self):
        return self*(-1.0)

    def __div__(self, other):
        try:
            return Point(self.x()/other, self.y()/other, self.z()/other)
        except ZeroDivisionError:
            print '** Error (division by 0 in Point.__div__).'
            return Point()

    def __str__(self):
        return '(%.2f, %.2f, %.2f)' % (self.x(), self.y(), self.z())



class Vector(Point):

    pass



def vector2point(v):
    return Point(v.x(), v.y(), v.z())



class CalMomentsData:

    def __init__(self, xtalData):
        self.XtalData = xtalData
        self.Point = vector2point(self.__getPoint())
        self.DistToAxis = 0.
        self.CoordAlongAxis = 0.
        self.XZMarker = ROOT.TMarker(self.Point.x(), self.Point.z(), 25)
        self.YZMarker = ROOT.TMarker(self.Point.y(), self.Point.z(), 25)
        self.XYMarker = ROOT.TMarker(self.Point.x(), self.Point.y(), 25)

    def setMaxWeight(self, maxWeight):
        try:
            markerSize = 1.8*sqrt(self.getWeight()/maxWeight)
        except ValueError:
            markerSize = 0.1
        markerSize = max(markerSize, 0.15)
        self.XZMarker.SetMarkerSize(markerSize)
        self.YZMarker.SetMarkerSize(markerSize)
        self.XYMarker.SetMarkerSize(markerSize)

    def setColor(self, color, redraw = False):
        self.XZMarker.SetMarkerColor(color)
        self.YZMarker.SetMarkerColor(color)
        self.XYMarker.SetMarkerColor(color)
        if redraw:
            self.XZMarker.Draw()
            self.YZMarker.Draw()
            self.XYMarker.Draw()

    def getTower(self):
        return self.XtalData.getPackedId().getTower()

    def getLayer(self):
        return self.XtalData.getPackedId().getLayer()

    def getColumn(self):
        return self.XtalData.getPackedId().getColumn()

    def getDistToPoint(self, point):
        return (self.Point - point).Mag()

    def __getPoint(self):
        return self.XtalData.getPosition()

    def getPoint(self):
        return self.Point

    def getWeight(self):
        return self.XtalData.getEnergy()

    def getDistToAxis(self):
        return self.DistToAxis

    def calcDistToAxis(self, centroid, axis):
        diffVec   = self.Point - centroid
        crossProd = -axis.Cross(diffVec)
        self.DistToAxis = crossProd.Mag()
        return self.DistToAxis

    def calcCoordAlongAxis(self, centroid, axis):
        self.CoordAlongAxis = (self.Point - centroid).Dot(axis)
        return self.CoordAlongAxis

    def __cmp__(self, other):
        if self.DistToAxis > other.DistToAxis:
            return 1
        else:
            return -1

    def __str__(self):
        return 'xtal (%d, %d, %d) with %.2f MeV @ %s' %\
            (self.getTower(), self.getLayer(), self.getColumn(),
             self.getWeight(), self.Point)
        


CAL_LAYOUT = CalLayout()


class CalClusterDisplay:

    def __init__(self, cluster, xtalList, maxWeight):
        # Setup the xtals...
        self.MomentsDataList = []
        centroid = cluster.getMomParams().getCentroid()
        axis = cluster.getMomParams().getAxis()
        for xtalData in xtalList:
            momentsData = CalMomentsData(xtalData)
            momentsData.calcCoordAlongAxis(centroid, axis)
            dist = momentsData.calcDistToAxis(centroid, axis)
            momentsData.setMaxWeight(maxWeight)
            self.MomentsDataList.append(momentsData)
        #... then the output from the moments analysis...
        momParams = cluster.getMomParams()
        xc = momParams.getCentroid().x()
        yc = momParams.getCentroid().y()
        zc = momParams.getCentroid().z()
        dx = momParams.getAxis().x()
        dy = momParams.getAxis().y()
        dz = momParams.getAxis().z()
        if dz == 0.:
            dz = 0.0001
        x1 = xc - (zc - Z_MIN)*(dx/dz)
        z1 = Z_MIN
        x2 = xc + (Z_MAX - zc)*(dx/dz)
        z2 = Z_MAX
        self.XZMomCentrMarker = ROOT.TMarker(xc, zc, 20)
        self.XZMomDirection = ROOT.TLine(x1, z1, x2, z2)
        self.XZMomDirection.SetLineWidth(1)
        y1 = yc - (zc - Z_MIN)*(dy/dz)
        y2 = yc + (Z_MAX - zc)*(dy/dz)
        self.YZMomCentrMarker = ROOT.TMarker(yc, zc, 20)
        self.YZMomDirection = ROOT.TLine(y1, z1, y2, z2)
        self.YZMomDirection.SetLineWidth(1)
        # and eventually the output from the fit.
        fitParams = cluster.getFitParams()
        xc = fitParams.getCentroid().x()
        yc = fitParams.getCentroid().y()
        zc = fitParams.getCentroid().z()
        dx = fitParams.getAxis().x()
        dy = fitParams.getAxis().y()
        dz = fitParams.getAxis().z()
        if dz == 0.:
            dz = 0.0001
        x1 = xc - (zc - Z_MIN)*(dx/dz)
        z1 = Z_MIN
        x2 = xc + (Z_MAX - zc)*(dx/dz)
        z2 = Z_MAX
        self.XZFitCentrMarker = ROOT.TMarker(xc, zc, 24)
        self.XZFitDirection = ROOT.TLine(x1, z1, x2, z2)
        self.XZFitDirection.SetLineWidth(1)
        self.XZFitDirection.SetLineStyle(7)
        y1 = yc - (zc - Z_MIN)*(dy/dz)
        y2 = yc + (Z_MAX - zc)*(dy/dz)
        self.YZFitCentrMarker = ROOT.TMarker(yc, zc, 24)
        self.YZFitDirection = ROOT.TLine(y1, z1, y2, z2)
        self.YZFitDirection.SetLineWidth(1)
        self.YZFitDirection.SetLineStyle(7)

    def setColor(self, color):
        for momentsData in self.MomentsDataList:
            momentsData.setColor(color)
        self.XZMomCentrMarker.SetMarkerColor(color)
        self.XZMomDirection.SetLineColor(color)
        self.YZMomCentrMarker.SetMarkerColor(color)
        self.YZMomDirection.SetLineColor(color)
        self.XZFitCentrMarker.SetMarkerColor(color)
        self.XZFitDirection.SetLineColor(color)
        self.YZFitCentrMarker.SetMarkerColor(color)
        self.YZFitDirection.SetLineColor(color)

    def drawXZ(self):
        for momentsData in self.MomentsDataList:
            momentsData.XZMarker.Draw()
        self.XZMomCentrMarker.Draw()
        self.XZMomDirection.Draw()
        self.XZFitCentrMarker.Draw()
        self.XZFitDirection.Draw()

    def drawYZ(self):
        for momentsData in self.MomentsDataList:
            momentsData.YZMarker.Draw()
        self.YZMomCentrMarker.Draw()
        self.YZMomDirection.Draw()
        self.YZFitCentrMarker.Draw()
        self.YZFitDirection.Draw()
        
    def drawXY(self):
        for momentsData in self.MomentsDataList:
            momentsData.XYMarker.Draw()

    def draw(self, sideCanvas, topCanvas):
        sideCanvas.cd(1)
        self.drawXZ()
        sideCanvas.cd(2)
        self.drawYZ()
        sideCanvas.cd()
        sideCanvas.Update()
        topCanvas.cd()
        self.drawXY()
        topCanvas.Update()
        


class CalDisplay(ReconReader):

    def __init__(self, reconFilePath, meritFilePath = None,
                 relFilePath = None, cut = '1'):
        ReconReader.__init__(self, reconFilePath, meritFilePath,
                             relFilePath, cut)
        self.SideCanvas = getSideCanvas('CalDisplay', 'CAL display')
        self.TopCanvas = getTopCanvas('CalTopDisplay', 'CAL top display')
        self.ClusterDisplayList = []

    def drawCalLayout(self):
        self.SideCanvas.cd(1)
        CAL_LAYOUT.draw('xz')
        self.SideCanvas.cd(2)
        CAL_LAYOUT.draw('yz')
        self.SideCanvas.Update()
        self.TopCanvas.cd()
        CAL_LAYOUT.draw('xy')
        self.TopCanvas.Update()

    def drawMcDir(self, color = ROOT.kGray, style = 1):
        xc = self.getMeritVariable('McX0')
        yc = self.getMeritVariable('McY0')
        zc = self.getMeritVariable('McZ0')
        dx = self.getMeritVariable('McXDir')
        dy = self.getMeritVariable('McYDir')
        dz = self.getMeritVariable('McZDir')
        if (xc is None) or (yc is None) or (zc is None):
            print "No MC info..."
            return None
        if dz == 0.:
            dz = 0.0001
        x1 = xc - (zc - Z_MIN)*(dx/dz)
        z1 = Z_MIN
        x2 = xc + (Z_MAX - zc)*(dx/dz)
        z2 = Z_MAX
        self.XZMcDirection = ROOT.TLine(x1, z1, x2, z2)
        self.XZMcDirection.SetLineColor(color)
        self.XZMcDirection.SetLineWidth(1)
        self.XZMcDirection.SetLineStyle(style)
        self.SideCanvas.cd(1)
        self.XZMcDirection.Draw()
        y1 = yc - (zc - Z_MIN)*(dy/dz)
        y2 = yc + (Z_MAX - zc)*(dy/dz)
        self.YZMcDirection = ROOT.TLine(y1, z1, y2, z2)
        self.YZMcDirection.SetLineColor(color)
        self.YZMcDirection.SetLineWidth(1)
        self.YZMcDirection.SetLineStyle(style)
        self.SideCanvas.cd(2)
        self.YZMcDirection.Draw()

    def drawCluster(self, i = 0):
        self.drawCalLayout()
        self.ClusterDisplayList = []
        cluster = self.getCalCluster(i)
        cluster.Print()
        xtalList = self.getCalClusterXtalList(i)
        uberCluster = self.getCalUberCluster()
        maxEnergy = uberCluster.getMSTreeParams().getMaxXtalEnergy()
        self.drawMcDir()
        clusterDisplay = CalClusterDisplay(cluster, xtalList, maxEnergy)
        clusterDisplay.draw(self.SideCanvas, self.TopCanvas)
        self.ClusterDisplayList.append(clusterDisplay)

    def drawClusters(self):
        for (i, cluster) in enumerate(self.getCalClusterCol()):
            print 'Cluster %d:' % i
            print '-'*80
            cluster.getXtalsParams().Print()
            print '-'*20
            cluster.getMSTreeParams().Print()
            print '-'*20
            cluster.getMomParams().Print()
            print '-'*20
            cluster.getClassParams().Print()
            print '-'*80
            print
        self.drawCalLayout()
        self.ClusterDisplayList = []
        uberCluster = self.getCalUberCluster()
        uberXtals = self.getCalUberClusterXtalList()
        maxEnergy = uberCluster.getMSTreeParams().getMaxXtalEnergy()
        self.drawMcDir(ROOT.kBlack)
        display = CalClusterDisplay(uberCluster, uberXtals, maxEnergy)
        if self.RelTable != None:
            display.setColor(ROOT.kGray)
        display.draw(self.SideCanvas, self.TopCanvas)
        self.ClusterDisplayList.append(display)
        
        nClusters = self.getNumClusters()
        if nClusters<=1:
            maxCluster = 1
        elif nClusters>len(COLOR_WHEEL):
            maxCluster = len(COLOR_WHEEL)
        else:
            maxCluster = nClusters-1        
        for (i, cluster) in enumerate(self.getCalClusterCol()):
            if i >= maxCluster:
                break
            xtals = self.getCalClusterXtalList(i)
            display = CalClusterDisplay(cluster, xtals, maxEnergy)
            display.setColor(COLOR_WHEEL[i])
            display.draw(self.SideCanvas, self.TopCanvas)
            self.ClusterDisplayList.append(display)




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
    display = CalDisplay(reconFilePath, meritFilePath, relFilePath, opts.c)
    answer = ''
    eventNumber = 0
    while answer != 'q':
        print 'ReconReader retrieving event %d...' % eventNumber
        cutPassed = display.getEntry(eventNumber)
        numClusters = display.getNumClusters()
        numXtals = display.getCalTotalNumXtals()
        print '%d cluster(s) found (including the uber), %d xtals in total.' %\
            (numClusters, numXtals)
        if numClusters > 0 and cutPassed:
            display.drawClusters()
            answer = raw_input('\n\n\tPress q to quit, s to save or type a number...')
        try:
            eventNumber = int(answer)
        except:
            eventNumber += 1
