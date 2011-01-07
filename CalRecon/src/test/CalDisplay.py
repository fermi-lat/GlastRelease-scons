
import math

from ReconReader import *


CSI_RAD_LEN = 18.6
CSI_MOL_RAD = 35.7

CAL_TOWER_GAP = TOWER_PITCH - 0.5*(CSI_WIDTH*12 + CSI_LENGTH)


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
        if markerSize > 0.15:
            self.XZMarker.SetMarkerSize(markerSize)
            self.YZMarker.SetMarkerSize(markerSize)
            self.XYMarker.SetMarkerSize(markerSize)
        else:
            self.XZMarker.SetMarkerStyle(1)
            self.YZMarker.SetMarkerStyle(1)
            self.XYMarker.SetMarkerStyle(1)

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

class CalDisplay(ReconReader):

    def __init__(self, reconFilePath, meritFilePath = None,
                 relFilePath = None, cut = '1'):
        ReconReader.__init__(self, reconFilePath, meritFilePath,
                             relFilePath, cut)
        self.Canvas = getFullCanvas('CalDisplay', 'CAL display')
        self.TopCanvas = getTopCanvas('CalTopDisplay', 'CAL top display')
        self.Cluster = None
        self.XtalEneHist = None

    def __refineLongProfile(self, transScale = 0.9):
        self.LongProfileRef = ROOT.TGraph()
        self.LongProfileRef.SetMarkerStyle(24)
        self.LongProfileRef.SetMarkerSize(0.8)
        axis = self.Cluster.getMomParams().getAxis()
        xdir = axis.x()
        ydir = axis.y()
        zdir = axis.z()
        weightSum = 0.
        tsum = 0.
        tsum2 = 0.
        tsum3 = 0.
        transRms = self.Cluster.getMomParams().getTransRms()
        n = 0
        tmin = 10000.
        tmax = -10000.
        for momentsData in self.MomentsDataList:
            dist = momentsData.DistToAxis
            if dist < transRms*transScale:
                t = momentsData.CoordAlongAxis
                # Get the tower Id and the position into the grid.
                towerId = momentsData.getTower()
                numGapsX = towerId % 4
                numGapsY = towerId / 4
                # Correct for the gaps (not exact!!)
                if abs(xdir) > 0.05:
                    t -= numGapsX*CAL_TOWER_GAP/xdir
                if abs(ydir) > 0.05:
                    t -= numGapsY*CAL_TOWER_GAP/ydir
                # Normalize to the X0 in CsI.
                t /= CSI_RAD_LEN
                if t < tmin:
                    tmin = t
                if t > tmax:
                    tmax = t
                energy = momentsData.getWeight()
                weightSum += energy
                tsum  += energy*t
                tsum2 += energy*t*t
                tsum3 += energy*t*t*t
                self.LongProfileRef.SetPoint(n, t, energy)
                n += 1
        if weightSum > 0:
            tsum  /= weightSum
            tsum2 /= weightSum
            tsum3 /= weightSum
            if tsum2 > 0.:
                mean = tsum
                var = tsum2 - mean*mean
                skew = tsum3 - 3*var*mean - mean*mean*mean
                if var > 0:
                    skew /= (var**1.5)
            
        self.Canvas.cd(2).cd(3)
        self.LongProfileRef.GetXaxis().SetTitle('Longitudinal position (X0)')
        self.LongProfileRef.GetXaxis().SetTitleOffset(2.00)
        self.LongProfileRef.GetYaxis().SetTitle('Xtal energy (MeV)')
        self.LongProfileRef.GetYaxis().SetTitleOffset(3.50)
        self.LongProfileRef.Draw('ap')

    def __drawXtals(self):
        maxWeight = self.Cluster.getMSTreeParams().getMaxXtalEnergy()
        centroid = self.Cluster.getMomParams().getCentroid()
        axis = self.Cluster.getMomParams().getAxis()
        self.MomentsDataList = []
        for xtalData in self.XtalList:
            momentsData = CalMomentsData(xtalData)
            momentsData.calcCoordAlongAxis(centroid, axis)
            dist = momentsData.calcDistToAxis(centroid, axis)
            momentsData.setMaxWeight(maxWeight)
            self.MomentsDataList.append(momentsData)
        # Draw xtals in XZ view.
        self.Canvas.cd(1).cd(1)
        CAL_LAYOUT.draw('xz')
        for momentsData in self.MomentsDataList:
            momentsData.XZMarker.Draw()
        # Draw xtals in YZ view.
        self.Canvas.cd(1).cd(2)
        CAL_LAYOUT.draw('yz')
        for momentsData in self.MomentsDataList:
            momentsData.YZMarker.Draw()
        # Draw xtals in the XY view.
        self.TopCanvas.cd()
        CAL_LAYOUT.draw('xy')
        for momentsData in self.MomentsDataList:
            momentsData.XYMarker.Draw()
        self.TopCanvas.Update()

    def __drawMcDir(self):
        xc = self.getMeritVariable('McX0')
        yc = self.getMeritVariable('McY0')
        zc = self.getMeritVariable('McZ0')
        dx = self.getMeritVariable('McXDir')
        dy = self.getMeritVariable('McYDir')
        dz = self.getMeritVariable('McZDir')
        if dz == 0.:
            dz = 0.0001
        x1 = xc - (zc - Z_MIN)*(dx/dz)
        z1 = Z_MIN
        x2 = xc + (Z_MAX - zc)*(dx/dz)
        z2 = Z_MAX
        self.XZMcDirection = ROOT.TLine(x1, z1, x2, z2)
        self.XZMcDirection.SetLineColor(ROOT.kBlue)
        self.XZMcDirection.SetLineWidth(1)
        self.XZMcDirection.SetLineStyle(7)
        self.Canvas.cd(1).cd(1)
        self.XZMcDirection.Draw()
        y1 = yc - (zc - Z_MIN)*(dy/dz)
        y2 = yc + (Z_MAX - zc)*(dy/dz)
        self.YZMcDirection = ROOT.TLine(y1, z1, y2, z2)
        self.YZMcDirection.SetLineColor(ROOT.kBlue)
        self.YZMcDirection.SetLineWidth(1)
        self.YZMcDirection.SetLineStyle(7)
        self.Canvas.cd(1).cd(2)
        self.YZMcDirection.Draw()

    def __drawCalMomParams(self):
        momParams = self.Cluster.getMomParams()
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
        self.XZCentrMarker = ROOT.TMarker(xc, zc, 20)
        self.XZCentrMarker.SetMarkerColor(ROOT.kRed)
        self.XZDirection = ROOT.TLine(x1, z1, x2, z2)
        self.XZDirection.SetLineColor(ROOT.kRed)
        self.XZDirection.SetLineWidth(1)
        self.Canvas.cd(1).cd(1)
        self.XZCentrMarker.Draw()
        self.XZDirection.Draw()
        y1 = yc - (zc - Z_MIN)*(dy/dz)
        y2 = yc + (Z_MAX - zc)*(dy/dz)
        self.YZCentrMarker = ROOT.TMarker(yc, zc, 20)
        self.YZCentrMarker.SetMarkerColor(ROOT.kRed)
        self.YZDirection = ROOT.TLine(y1, z1, y2, z2)
        self.YZDirection.SetLineColor(ROOT.kRed)
        self.YZDirection.SetLineWidth(1)
        self.Canvas.cd(1).cd(2)
        self.YZCentrMarker.Draw()
        self.YZDirection.Draw()

    def __drawXtalEneDist(self, truncFact = 0.02):
        if self.XtalEneHist is not None:
            self.XtalEneHist.Delete()
        maxWeight = self.Cluster.getMSTreeParams().getMaxXtalEnergy()
        totalEnergy = self.Cluster.getXtalsParams().getXtalRawEneSum()
        self.XtalEneHist = ROOT.TH1F('xtalEneDist', 'Xtal energy distribution',
                                     25, 0, 1.1*maxWeight)
        for xtalData in self.XtalList:
            xtalEnergy = xtalData.getEnergy()
            if xtalEnergy > truncFact*totalEnergy:
                self.XtalEneHist.Fill(xtalEnergy)
        self.Canvas.cd(2).cd(1)
        ROOT.gPad.SetLogy(True)
        self.XtalEneHist.SetXTitle('Xtal energy (MeV)')
        self.XtalEneHist.GetXaxis().SetTitleOffset(2.00)
        self.XtalEneHist.Draw()

    def __drawLongProfile(self, transScale = 0.9):
        self.LongProfile = ROOT.TGraph()
        self.LongProfile.SetMarkerStyle(24)
        self.LongProfile.SetMarkerSize(0.8)
        transRms = self.Cluster.getMomParams().getTransRms()
        i = 0
        for momentsData in self.MomentsDataList:
            dist = momentsData.DistToAxis
            if dist < transRms*transScale:
                self.LongProfile.SetPoint(i, momentsData.CoordAlongAxis,
                                          momentsData.getWeight())
                i += 1
        self.Canvas.cd(2).cd(2)
        self.LongProfile.GetXaxis().SetTitle('Longitudinal position (mm)')
        self.LongProfile.GetXaxis().SetTitleOffset(2.00)
        self.LongProfile.GetYaxis().SetTitle('Xtal energy (MeV)')
        self.LongProfile.GetYaxis().SetTitleOffset(3.50)
        self.LongProfile.Draw('ap')

    def drawCluster(self, i = 0):
        print 'Full info for the first cluster...'
        self.Cluster = self.getCalCluster(i)
        self.XtalList = self.getCalClusterXtalList(i)
        self.Cluster.Print()
        # Draw the CAL layout and the actual xtals.
        self.__drawXtals()
        # Draw direction/centroid from the moments analysis.
        self.__drawCalMomParams()
        # Draw the MC directions.
        self.__drawMcDir()
        # Draw the various histograms/graphs.
        self.__drawXtalEneDist()
        self.__drawLongProfile()
        self.__refineLongProfile()
        # Update canvas.
        self.Canvas.cd()
        self.Canvas.Update()



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
        display.getEntry(eventNumber)
        numClusters = display.getNumClusters()
        numXtals = display.getCalTotalNumXtals()
        print '%d cluster(s) found, %d xtals in total.' %\
            (numClusters, numXtals)
        if numClusters > 0:
            display.drawCluster(0)
            answer = raw_input('Press q to quit, s to save or type a number...')
        try:
            eventNumber = int(answer)
        except:
            eventNumber += 1
