
from ReconReader import *


class Point(ROOT.TVector3):

    def __init__(self, x = 0., y = 0., z = 0.):
        ROOT.TVector3.__init__(self, x, y, z)

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
        return '(%.4f, %.4f, %.4f)' % (self.x(), self.y(), self.z())



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

    def setMaxWeight(self, maxWeight):
        markerSize = 1.8*sqrt(self.getWeight()/maxWeight)
        if markerSize > 0.15:
            self.XZMarker.SetMarkerSize(markerSize)
            self.YZMarker.SetMarkerSize(markerSize)
        else:
            self.XZMarker.SetMarkerStyle(1)
            self.YZMarker.SetMarkerStyle(1)

    def setX(self, x):
        y = self.Point.y()
        z = self.Point.z()
        self.Point = Point(x, y, z)
        self.XZMarker = ROOT.TMarker(self.Point.x(), self.Point.z(), 25)

    def setY(self, y):
        x = self.Point.x()
        z = self.Point.z()
        self.Point = Point(x, y, z)
        self.YZMarker = ROOT.TMarker(self.Point.y(), self.Point.z(), 25) 

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
        diffVec   = centroid - self.Point
        crossProd = axis.Cross(diffVec)
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


CAL_LAYOUT = CalLayout()

class CalDisplay(ReconReader):

    def __init__(self, reconFilePath, meritFilePath = None,
                 relFilePath = None, cut = '1'):
        ReconReader.__init__(self, reconFilePath, meritFilePath,
                             relFilePath, cut)
        self.Canvas = getSideCanvas('CalDisplay', 'CAL display')

    def drawCluster(self, i = 0):
        print 'Full info for the first cluster...'
        cluster = self.getCalCluster(i)
        cluster.Print()
        maxWeight = cluster.getMSTreeParams().getMaxXtalEnergy()
        xtalList = self.getCalClusterXtalList(i)
        self.MomentsDataList = []
        for xtalData in xtalList:
            momentsData = CalMomentsData(xtalData)
            momentsData.setMaxWeight(maxWeight)
            self.MomentsDataList.append(momentsData)
        # Draw XZ view.
        self.Canvas.cd(1)
        CAL_LAYOUT.draw('xz')
        for momentsData in self.MomentsDataList:
            momentsData.XZMarker.Draw()

        # Draw YZ view.
        self.Canvas.cd(2)
        CAL_LAYOUT.draw('yz')
        for momentsData in self.MomentsDataList:
            momentsData.YZMarker.Draw()

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
