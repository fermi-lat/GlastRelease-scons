
# Author(s):
# Luca Baldini (luca.baldini@pi.infn.it)
# Johan Bregeon (johan.bregeon@pi.infn.it)
#
# Purpose:
# This is essentially a python wrapper of the CAL moments analysy which
# is potentially useful when experimenting with the algorithms (you don't have
# to recompile and run GR).
# It has the capability to play-back recon files and automatically compares
# the output of the python code with the numbers written in the file itself
# at creation time.


import os
import sys
import ROOT
import copy
import numpy

from math import acos, sqrt, cos, pi, fabs, log

from CalLayout   import *
from ReconReader import *


ROOT.gStyle.SetCanvasColor(ROOT.kWhite)

M_PI = pi
MOLIERE_RADIUS = 35.3

LIBRARIES = ['libcommonRootData.so', 'libreconRootData.so']
MERIT_VARS = ['EvtEventId',
              'McEnergy', 'CalEnergyRaw', 'CTBBestEnergy', 'CalEnergyCorr',
              'McXDir', 'McYDir', 'McZDir',
              'Tkr1XDir', 'Tkr1YDir', 'Tkr1ZDir',
              'Tkr1X0', 'Tkr1Y0', 'Tkr1Z0',
              'CalCsIRLn', 'CalLATRLn',
              'CalTrackDoca', 'CalTrackAngle', 'CTBCORE',
              'CalTransRms', 'CalLongRms', 'CalLRmsAsym',
              'TkrNumTracks']
MIN_NUM_XTALS = 3

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


CSI_CRITICAL_ENERGY = 11.17

def Gamma(a):
    return ROOT.TMath.Gamma(a)

def gamma(a, x):
    return ROOT.TMath.Gamma(a, x)

def GGamma(a, b, t1, t2):
    return (Gamma(a)/(b**a))*(gamma(a, b*t2) - gamma(a, b*t1))

def getExpectedLongParameters(energy, tmin, tmax, c = -0.5):
    if energy <= 0 or tmax <= tmin:
        return (-1, -1, -1)
    b     = 0.5
    a     = b*(log(energy/CSI_CRITICAL_ENERGY) + c)
    mean = GGamma(a+2, b, tmin, tmax)/GGamma(a+1, b, tmin, tmax)
    sigma2 = GGamma(a+3, b, tmin, tmax)/GGamma(a+1, b, tmin, tmax) -\
             mean**2
    skewness = (GGamma(a+4, b, tmin, tmax)/GGamma(a+1, b, tmin, tmax) -\
                3*mean*sigma2 - mean**3)/(sigma2**1.5)
    # Change sign to the skewness, as the axis points up!
    return (mean, sigma2, -skewness)


M_IDENTITY_3_3 = numpy.identity(3)
M_ZEROS_3_3    = numpy.zeros((3,3))
M_ZEROS_9_3    = numpy.zeros((9,3))



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


# Convenience class to store the information about the single iterations
# of the moments analysis without changing the too much the code in
# the CalMomentsAnalysis class.

CAL_LAYOUT = CalLayout()


class CalMomentsAnalysisIteration:

    def __init__(self, momentsAnalysis, dataVec):
        self.MomentsAnalysis = momentsAnalysis
        self.IterationNumber = copy.deepcopy(momentsAnalysis.NumIterations)
        self.NumXtals = len(dataVec)
        self.Centroid = copy.deepcopy(momentsAnalysis.Centroid)
        self.Centroid = vector2point(self.Centroid)
        self.Moment = copy.deepcopy(momentsAnalysis.Moment)
        self.RmsLong = copy.deepcopy(momentsAnalysis.RmsLong)
        self.RmsTrans = copy.deepcopy(momentsAnalysis.RmsTrans)
        self.RmsLongAsym = copy.deepcopy(momentsAnalysis.RmsLongAsym)
        self.WeightSum = copy.deepcopy(momentsAnalysis.WeightSum)
        self.SkewnessLong = copy.deepcopy(momentsAnalysis.SkewnessLong)
        # Normalize the skewness.
        sigma = sqrt(self.RmsLong/self.WeightSum)
        self.SkewnessLong = self.SkewnessLong/((sigma**3)*self.WeightSum)
        self.NumIterations = copy.deepcopy(momentsAnalysis.NumIterations)
        self.Axis = copy.deepcopy(momentsAnalysis.Axis)
        for (i, axis) in enumerate(self.Axis):
            self.Axis[i] = vector2point(axis)
        self.WeightSum = copy.deepcopy(momentsAnalysis.WeightSum)
        self.InertiaTensor = copy.deepcopy(momentsAnalysis.InertiaTensor)
        self.PrincipalAxis = copy.deepcopy(momentsAnalysis.PrincipalAxis)
        # Construct the longitudinal profile
        # (and store the CAL moments data)
        self.LongProfile = ROOT.TGraph()
        self.LongProfile.SetMarkerStyle(24)
        self.LongProfile.SetMarkerSize(0.8)
        self.MomentsDataList = []
        self.MomentsDataListPP = []
        for (i, dataPoint) in enumerate(dataVec):
            self.LongProfile.SetPoint(i, dataPoint.CoordAlongAxis,
                                      dataPoint.getWeight())
            dataPoint.setMaxWeight(momentsAnalysis.MaxWeight)
            # Need a copy/deepcopy in the following 2 lines?
            dp = copy.deepcopy(dataPoint)
            dp.Point = vector2point(dp.Point)
            self.MomentsDataList.append(dp)
            dp = copy.deepcopy(dataPoint)
            dp.Point = vector2point(dp.Point)
            self.MomentsDataListPP.append(dp)
        self.CalTkrAngle = -1
        self.calculateLayerCentroids()

    def calculateLayerCentroids(self):
        # Calculate the layer centroids.
        self.XZLayerCentrDict = {}
        self.YZLayerCentrDict = {}
        self.LayerWeightDict  = {}
        self.XZLayerCentrMarkerDict = {}
        self.YZLayerCentrMarkerDict = {}
        for layer in range(8):
            self.XZLayerCentrDict[layer] = 0.0
            self.YZLayerCentrDict[layer] = 0.0
            self.LayerWeightDict[layer]  = 0.0
        for momentsData in self.MomentsDataList:
            layer = momentsData.getLayer()
            weight = momentsData.getWeight()
            self.XZLayerCentrDict[layer] += momentsData.Point.x()*weight
            self.YZLayerCentrDict[layer] += momentsData.Point.y()*weight
            self.LayerWeightDict[layer]  += weight
        for layer in range(8):
            try:
                self.XZLayerCentrDict[layer] /= self.LayerWeightDict[layer]
                self.YZLayerCentrDict[layer] /= self.LayerWeightDict[layer]
                z = getCalLayerZ(layer)
                mStyle = 5
                mSize = 1
                mxz = ROOT.TMarker(self.XZLayerCentrDict[layer], z, mStyle)
                mxz.SetMarkerSize(mSize)
                myz = ROOT.TMarker(self.YZLayerCentrDict[layer], z, mStyle)
                myz.SetMarkerSize(mSize)
                self.XZLayerCentrMarkerDict[layer] = mxz
                self.YZLayerCentrMarkerDict[layer] = myz
            except:
                pass
        # Now a look at residuals &Co.
        resThreshold = 15.0 #mm
        fitFunc = ROOT.TF1('fitFunc', 'pol1')
        self.XResDict = {}
        self.YResDict = {}
        self.XPPLayerList = []
        self.YPPLayerList = []
        for layer in range(8):
            if self.LayerWeightDict[layer]:
                # First get the x and y positions that you would expect from
                # the other layers.
                gxz = ROOT.TGraphErrors()
                gyz = ROOT.TGraphErrors()
                nxz = 0
                for i in range(8):
                    w = self.LayerWeightDict[i]
                    if w > 0 and i != layer:
                        x = self.XZLayerCentrDict[i]
                        z = getCalLayerZ(i)
                        gxz.SetPoint(nxz, z, x)
                        gxz.SetPointError(nxz, 0, sqrt(w/self.WeightSum))
                        nxz += 1
                nyz = 0
                for i in range(8):
                    w = self.LayerWeightDict[i]
                    if w > 0 and i != layer:
                        y = self.YZLayerCentrDict[i]
                        z = getCalLayerZ(i)
                        gyz.SetPoint(nyz, z, y)
                        gyz.SetPointError(nyz, 0, sqrt(w/self.WeightSum))
                        nyz += 1
                x = self.XZLayerCentrDict[layer]
                y = self.YZLayerCentrDict[layer]
                z = getCalLayerZ(layer)
                gxz.Fit('fitFunc', 'Q')
                xfit = fitFunc.Eval(z)
                gyz.Fit('fitFunc', 'Q')
                yfit = fitFunc.Eval(z)
                self.XResDict[layer] = x - xfit
                self.YResDict[layer] = y - yfit
                if abs(self.XResDict[layer]) > resThreshold:
                    if layer in [0, 2, 4, 6]:
                        mxz = self.XZLayerCentrMarkerDict[layer]
                        mxz.SetMarkerColor(ROOT.kRed)
                        self.XPPLayerList.append(layer)
                if abs(self.YResDict[layer]) > resThreshold:
                    if layer in [1, 3, 5, 7]:
                        myz = self.YZLayerCentrMarkerDict[layer]
                        myz.SetMarkerColor(ROOT.kRed)
                        self.YPPLayerList.append(layer)
        # And eventually postprocess the dataVec.
        for momentsData in self.MomentsDataListPP:
            layer = momentsData.getLayer()
            if layer in self.XPPLayerList:
                x = momentsData.Point.x() - self.XResDict[layer]
                momentsData.setX(x)
            if layer in self.YPPLayerList:
                y = momentsData.Point.y() - self.YResDict[layer]
                momentsData.setY(y)

    def draw(self, reconReader, layerCentroids = True):
        dx = self.Axis[1].x()
        dy = self.Axis[1].y()
        dz = self.Axis[1].z()
        xc = self.Centroid.x()
        yc = self.Centroid.y()
        zc = self.Centroid.z()
        cluster = reconReader.getCalRecon().getCalClusterCol().At(0)
        rxc = cluster.getParams().getCentroid().x()
        ryc = cluster.getParams().getCentroid().y()
        rzc = cluster.getParams().getCentroid().z()
        rdx = cluster.getParams().getAxis().x()
        rdy = cluster.getParams().getAxis().y()
        rdz = cluster.getParams().getAxis().z()
        if reconReader.MeritChain is not None:
            tdx = reconReader.getMeritVariable('Tkr1XDir')
            tdy = reconReader.getMeritVariable('Tkr1YDir')
            tdz = reconReader.getMeritVariable('Tkr1ZDir')
            txc = reconReader.getMeritVariable('Tkr1X0')
            tyc = reconReader.getMeritVariable('Tkr1Y0')
            tzc = reconReader.getMeritVariable('Tkr1Z0')
            tkrDir = Vector(tdx, tdy, tdz)
            calDir = self.Axis[1]
            self.CalTkrAngle = (180./M_PI)*calDir.Angle(-tkrDir)
        cName = 'cMomIter%d' % self.IterationNumber
        cTitle = 'Moments analysis---iteration %d' % self.IterationNumber
        self.Canvas = getSideCanvas(cName, cTitle)
        self.Canvas.cd(1)
        # Draw XZ view.
        CAL_LAYOUT.draw('xz')
        text = 'CAL-TKR 3d angle: %.2f#circ (original %.2f#circ)' %\
               (self.CalTkrAngle,
                (180./M_PI)*reconReader.getMeritVariable('CalTrackAngle'))
        self.CalTkrAngleLabel = ROOT.TLatex(0.7, 0.92, text)
        self.CalTkrAngleLabel.SetNDC()
        self.CalTkrAngleLabel.SetTextSize(0.08)
        self.CalTkrAngleLabel.Draw()
        # Draw xtals.
        for dataPoint in  self.MomentsDataList:
            dataPoint.XZMarker.Draw()
        # Draw layer centroids.
        if layerCentroids:
            for marker in self.XZLayerCentrMarkerDict.values():
                marker.Draw()
        # Draw centroid from the python code.
        self.XZCentrMarker = ROOT.TMarker(xc, zc, 20)
        self.XZCentrMarker.SetMarkerColor(ROOT.kRed)
        self.XZCentrMarker.Draw()
        # Draw direction from the python code.
        x1 = xc - (zc - Z_MIN)*(dx/dz)
        z1 = Z_MIN
        x2 = xc + (Z_MAX - zc)*(dx/dz)
        z2 = Z_MAX
        self.XZDirection = ROOT.TLine(x1, z1, x2, z2)
        self.XZDirection.SetLineColor(ROOT.kRed)
        self.XZDirection.SetLineWidth(1)
        self.XZDirection.Draw()
        # If the numbers are different wrt the recon file, draw the
        # recon information.
        if not self.MomentsAnalysis.EqualToRecon:
            # Draw centroid from the recon file.
            self.XZReconCentrMarker = ROOT.TMarker(rxc, rzc, 30)
            self.XZReconCentrMarker.SetMarkerColor(ROOT.kRed)
            self.XZReconCentrMarker.Draw()
            # Draw direction from the recon file.
            rx1 = rxc - (rzc - Z_MIN)*(rdx/rdz)
            rz1 = Z_MIN
            rx2 = rxc + (Z_MAX - rzc)*(rdx/rdz)
            rz2 = Z_MAX
            self.XZReconDirection = ROOT.TLine(rx1, rz1, rx2, rz2)
            self.XZReconDirection.SetLineColor(ROOT.kRed)
            self.XZReconDirection.SetLineStyle(7)
            self.XZReconDirection.SetLineWidth(1)
            self.XZReconDirection.Draw()
        # If we do have the Merit, draw the TKR direction
        if reconReader.MeritChain is not None:
            tx1 = txc - (tzc - Z_MIN)*(tdx/tdz)
            tz1 = Z_MIN
            tx2 = txc + (Z_MAX - tzc)*(tdx/tdz)
            tz2 = Z_MAX
            self.XZTkrDirection = ROOT.TLine(tx1, tz1, tx2, tz2)
            self.XZTkrDirection.SetLineColor(ROOT.kBlue)
            self.XZTkrDirection.SetLineWidth(1)
            self.XZTkrDirection.SetLineStyle(7)
            self.XZTkrDirection.Draw()
        self.Canvas.cd(2)
        # Draw YZ view.
        CAL_LAYOUT.draw('yz')
        # Draw xtals.
        for dataPoint in  self.MomentsDataList:
            dataPoint.YZMarker.Draw()
        # Draw layer centroids.
        if layerCentroids:
            for marker in self.YZLayerCentrMarkerDict.values():
                marker.Draw()
        # Draw centroid from the python code.
        self.YZCentrMarker = ROOT.TMarker(yc, zc, 20)
        self.YZCentrMarker.SetMarkerColor(ROOT.kRed)
        self.YZCentrMarker.Draw()
        # Draw direction from the python code.
        y1 = yc - (zc - Z_MIN)*(dy/dz)
        y2 = yc + (Z_MAX - zc)*(dy/dz)
        self.YZDirection = ROOT.TLine(y1, z1, y2, z2)
        self.YZDirection.SetLineColor(ROOT.kRed)
        self.YZDirection.SetLineWidth(1)
        self.YZDirection.Draw()
        # If the numbers are different wrt the recon file, draw the
        # recon information.
        if not self.MomentsAnalysis.EqualToRecon:
            # Draw centroid from the recon file.
            self.YZReconCentrMarker = ROOT.TMarker(ryc, rzc, 30)
            self.YZReconCentrMarker.SetMarkerColor(ROOT.kRed)
            self.YZReconCentrMarker.Draw()
            # Draw direction from the recon file.
            ry1 = ryc - (rzc - Z_MIN)*(rdy/rdz)
            ry2 = ryc + (Z_MAX - rzc)*(rdy/rdz)
            self.YZReconDirection = ROOT.TLine(ry1, rz1, ry2, rz2)
            self.YZReconDirection.SetLineColor(ROOT.kRed)
            self.YZReconDirection.SetLineStyle(7)
            self.YZReconDirection.SetLineWidth(1)
            self.YZReconDirection.Draw()
        if reconReader.MeritChain is not None:
            ty1 = tyc - (tzc - Z_MIN)*(tdy/tdz)
            ty2 = tyc + (Z_MAX - tzc)*(tdy/tdz)
            self.YZTkrDirection = ROOT.TLine(ty1, tz1, ty2, tz2)
            self.YZTkrDirection.SetLineColor(ROOT.kBlue)
            self.YZTkrDirection.SetLineWidth(1)
            self.YZTkrDirection.SetLineStyle(7)
            self.YZTkrDirection.Draw()
        self.Canvas.cd()
        self.Canvas.Update()

    def drawLongProfile(self, expSkewness = None):
        self.LongProfile.Draw('ap')
        self.LongProfile.GetXaxis().SetTitle('Position along princ. axis (mm)')
        self.LongProfile.GetYaxis().SetTitle('Energy (MeV)')
        self.LongProfile.GetYaxis().SetTitleOffset(1.4)
        self.TextBox = ROOT.TPaveText(0.6, 0.75, 0.99, 0.99, 'NDC')
        self.TextBox.SetTextAlign(12)
        self.TextBox.AddText('Iteration %d (%d xtals)' % (self.IterationNumber,
                                                          self.NumXtals))
        self.TextBox.AddText('Skewness = %.3f' % self.SkewnessLong)
        if expSkewness is not None:
            normSkewness = self.SkewnessLong/expSkewness
            self.TextBox.AddText('Exp. Skewness = %.3f' % expSkewness)
            self.TextBox.AddText('Norm. skewness = %.3f' % normSkewness)
        self.TextBox.Draw()
        ROOT.gPad.Update()

    def save(self, filePath):
        epsFilePath = '%s.eps' % filePath
        pngFilePath = '%s.png' % filePath
        try:
            self.Canvas.SaveAs(epsFilePath)
            os.system('epstopdf %s' % epsFilePath)
            self.Canvas.SaveAs(pngFilePath)
        except:
            pass

    def __str__(self):
        text = '** Iteration %d (%d xtals) **\n' % (self.IterationNumber,
                                                    self.NumXtals)
        text += 'Centroid         = %s\n' % self.Centroid
        text += '-----------------------------------------------------------\n'
        text += 'Inertia tensor: \n%s\n' % self.InertiaTensor
        text += '-----------------------------------------------------------\n'
        text += 'Moments          = %s\n' % self.Moment
        for i in range(3):
            text += 'Axis %d           = %s\n' % (i, self.Axis[i])
        text += 'Longitudinal RMS = %.3e\n' % self.RmsLong
        text += 'Transverse RMS   = %.3e\n' % self.RmsTrans
        text += 'Long. asymmetry  = %f\n' % self.RmsLongAsym
        text += 'Long skewness    = %.3e\n' % self.SkewnessLong
        text += 'CalTkrAngle      = %.3f degrees' % self.CalTkrAngle
        return text

    

class CalMomentsAnalysis:

    def __init__(self, parent, clipDataVec):
        self.Parent = parent
        self.ClipDataVec = clipDataVec
        self.Centroid = Point()
        self.Moment = [0., 0., 0.]
        self.RmsLong = 0.
        self.RmsTrans = 0.
        self.RmsLongAsym = 0.
        self.SkewnessLong = 0.0
        self.NumIterations = 0
        self.NumDroppedPoints = 0
        self.Axis = [Vector(0., 0., 0.),
                     Vector(0., 0., 1.),
                     Vector(0., 0., 0.)]        
        self.WeightSum = 0.0
        self.MaxWeight = 0.0
        # New variables, for debugging purposes.
        self.InertiaTensor = None
        self.PrincipalAxis = None
        self.IterationList = []

    def getNumIterations(self):
        return len(self.IterationList)

    def getMomentsCentroid(self):
        return self.Centroid

    def getMomentsAxis(self):
        return self.Axis[1]

    def getLongitudinalRms(self):
        return self.RmsLong

    def getTransverseRms(self):
        return self.RmsTrans

    def getLongAsymmetry(self):
        return self.RmsLongAsym
    
    def getLongSkewness(self):
        return self.SkewnessLong

    def getNumIterations(self):
        return self.NumIterations

    def getNumDroppedPoints(self):
        return self.NumDroppedPoints

    def estimateCentroid(self, dataVec):
        print 'Estimating centroid...'
        weightSum = 0.
        centroid = Point()
        for dataPoint in dataVec:
            weight = dataPoint.getWeight()
            weightSum += weight
            centroid += vector2point(dataPoint.getPoint()) * weight
        centroid /= weightSum
        return centroid

    def doMomentsAnalysis(self, dataVec, iniCentroid):
        print 'Starting moments analysis (iteration %d, %d xtals)...' %\
              (self.NumIterations, len(dataVec))
        self.WeightSum = 0.
        chiSquare = -1.
        if len(dataVec) < 2:
            print 'Less than 2 xtals found, returning chi square = -1.'
            return chiSquare
        Ixx = 0.0  
        Iyy = 0.0  
        Izz = 0.0
        Ixy = 0.0  
        Ixz = 0.0  
        Iyz = 0.0
        weightSum = 0.0
        skewness = -9999.0
        centroid = Point(0., 0., 0.)
        for dataPoint in dataVec:
            # Construct elements of (symmetric) "Inertia" Tensor:
            # See Goldstein, 1965, Chapter 5 (especially, eqs. 5-6,7,22,26).
            # Analysis easy when translated to energy centroid.
            #  get pointer to the reconstructed data for given crystal
            weight = dataPoint.getWeight()
            hit = vector2point(dataPoint.getPoint()) - iniCentroid
            Rsq  = hit.mag2()
            xprm = hit.x()
            yprm = hit.y()
            zprm = hit.z()
            Ixx += (Rsq - xprm*xprm) * weight
            Iyy += (Rsq - yprm*yprm) * weight
            Izz += (Rsq - zprm*zprm) * weight
            Ixy -= xprm*yprm * weight
            Ixz -= xprm*zprm * weight
            Iyz -= yprm*zprm * weight
            weightSum += weight
            if weight > self.MaxWeight:
                self.MaxWeight = weight
            centroid  += vector2point(dataPoint.getPoint()) * weight
        # Render determinant of Inertia Tensor into cubic form.
        p = - (Ixx + Iyy + Izz)
        q =   Iyy*Izz + Iyy*Ixx + Izz*Ixx - (Ixy*Ixy + Iyz*Iyz + Ixz*Ixz)
        r = - Ixx*Iyy*Izz + Ixx*Iyz*Iyz + Iyy*Ixz*Ixz +\
            Izz*Ixy*Ixy - 2.*Ixy*Iyz*Ixz
        # See CRC's Standard Mathematical Tables (19th edition), pp 103-105.
        # The substitution, y = x - p/3 converts  y^3 + p*y^2 + q*y + r = 0
        # to the form  x^3 + a*x + b = 0 .  Then, if b^2/4 + a^3/27 < 0 ,
        # there will be three real roots -- guaranteed since the Inertia Tensor
        # is symmetric.  A second substitution, x = m*cos(psi) ,
        # yields the roots.
        a = (3.*q - p*p)/3.
        b = (2.*p*p*p - 9.*p*q + 27.*r)/27.
        rad_test = b*b/4. + a*a*a/27.
        if (rad_test < 0.) and  (Ixy != 0.) and (Ixz != 0.) and (Iyz != 0.):
            # Update the weight and centroid
            self.WeightSum  = weightSum
            self.Centroid   = centroid
            self.Centroid  /= weightSum
            # Construct the roots, which are the principal moments.
            m = 2. * sqrt(-a/3.)
            psi = acos( 3.*b/(a*m) ) / 3.
            self.Moment[0] = m * cos(psi) - p/3.
            self.Moment[1] = m * cos(psi + 2.*M_PI/3.) - p/3.
            self.Moment[2] = m * cos(psi + 4.*M_PI/3.) - p/3.
            # Construct direction cosines; dircos for middle root is parallel
            # to longest principal axis.
            for iroot in range(3): 
                A = Iyz * (Ixx - self.Moment[iroot]) - Ixy*Ixz
                B = Ixz * (Iyy - self.Moment[iroot]) - Ixy*Iyz
                C = Ixy * (Izz - self.Moment[iroot]) - Ixz*Iyz
                D = sqrt( 1. / ( 1./(A*A) + 1./(B*B) + 1./(C*C) ) ) / C
                self.Axis[iroot] = Vector(D*C/A, D*C/B, D)
                # Set axis to "point up"
                if self.Axis[iroot].z() < 0.:
                    self.Axis[iroot] = -self.Axis[iroot]
            # "Chi-squared" = sum of residuals about principal axis,
            # through centroid, using input weight
            chiSquare = 0.
            # Skewness = third moment along the main axis w.r.t. the energy
            # centroid.
            skewness = 0.
            for dataPoint in dataVec:
                distToAxis = dataPoint.calcDistToAxis(self.Centroid,
                                                      self.Axis[1])
                coordAlongAxis = dataPoint.calcCoordAlongAxis(self.Centroid,
                                                              self.Axis[1])
                chiSquare += dataPoint.getWeight() * (distToAxis**2.0)
                skewness += dataPoint.getWeight() * (coordAlongAxis**3.0)
            # Scale by number of data points
            chiSquare /= weightSum * len(dataVec)
            # Final calculations to return moment of principal axis and
            # average of other two
            longMag1 = fabs(self.Moment[0])
            longMag2 = fabs(self.Moment[2]) 
            self.RmsLong = (longMag1 + longMag2) / 2.
            self.RmsTrans =  fabs(self.Moment[1])
            self.RmsLongAsym = (longMag1 - longMag2)/(longMag1 + longMag2)
            self.SkewnessLong = skewness
            # Fill additional variables for debugging purposes.
            self.InertiaTensor = numpy.matrix([[Ixx, Ixy, Ixz],
                                               [Ixy, Iyy, Iyz],
                                               [Ixz, Iyz, Izz]],
                                              dtype = 'd')
            self.PrincipalAxis = numpy.matrix([[self.Axis[1].X()],
                                               [self.Axis[1].Y()],
                                               [self.Axis[1].Z()]],
                                              dtype = 'd') 
            self.IterationList.append(CalMomentsAnalysisIteration(self,
                                                                  dataVec))
        else:
            chiSquare = -1.
        print 'Done, returning chi square = %.3f.' % chiSquare
        return chiSquare

    def doIterativeMomentsAnalysis(self, dataVec, inputCentroid = None,
                                   scaleFactor = 1.0):
        chiSq = -1.
        iterate = True
        self.NumIterations = 0
        self.NumDroppedPoints = 0
        # Set the data centroid
        centroid = inputCentroid or self.estimateCentroid(dataVec)
        # Iterate until either failure (chiSq < 0) or all remaining data
        # points are within "range" of axis
        while iterate:
            # Do the standard moments analysis
            localChiSq = self.doMomentsAnalysis(dataVec, centroid);
            # Make sure it didn't fail on this iteration
            if localChiSq < 0.: 
                break
            # Update global chi-square to pick up this iteration's value
            chiSq = localChiSq
            # Update the centroid for subsequent passes
            centroid = self.getMomentsCentroid()
            # Get the transverse moment
            rmsTrans = self.getTransverseRms()
            # Convert to distance by scaling with the sum of the weights 
            # and taking the square root
            rmsTrans = sqrt(rmsTrans / self.WeightSum)
            # Sort the data points by their distance from the principal axis
            dataVec.sort()
            # Assume all data within range
            iterate = False
            self.NumIterations += 1
            # Check the last element in the now sorted list of data points
            # and see if it is out of range
            while (len(dataVec)):
                momentsData = dataVec[-1]
                # If out of range drop this point and loop back to check again
                if (momentsData.getDistToAxis() > scaleFactor * rmsTrans):
                    dataVec.pop()
                    iterate = True
                    self.NumDroppedPoints += 1    
                else:
                    break
            # Make it harder to drop points on each iteration by boosting the
            # scale factor
            scaleFactor *= 2.
            # Finally, make sure have enough points remaining to proceed!
            if (len(dataVec) < 3):
                break
            if self.ClipDataVec:
                dataVec = self.Parent.getClippedDataVec(dataVec)
        return (chiSq, dataVec)

    def getCovarianceMatrix(self, dataVec):
        (L0, L1, L2) = self.Moment
        #print 'Moments: %s' % self.Moment
        #print 'Inertia tensor:\n%s' % self.InertiaTensor
        #print self.Moment
	S = numpy.matrix([ [self.Axis[j][i] for i in range(3)] \
                           for j in range(3) ], dtype = 'd')
        #print self.InertiaTensor
        #print S * self.InertiaTensor * S.I
        #raw_input()
	# Starting to define things
        #print 'S =\n%s' % S
        S_p = numpy.matrix([ [ 0     ,  S[2,0], -S[1,0]],
                             [-S[2,0],  0     ,  S[0,0]],
                             [ S[1,0], -S[0,0],  0     ],
                             [ 0     ,  S[2,1], -S[1,1]],
                             [-S[2,1],  0     ,  S[0,1]],
                             [ S[1,1], -S[0,1],  0     ],
                             [ 0     ,  S[2,2], -S[1,2]],
                             [-S[2,2],  0     ,  S[0,2]],
                             [ S[1,2], -S[0,2],  0     ] ], dtype = 'd')

        #print 'Sp =\n%s' % S_p
        D = numpy.matrix([ [1, 0,   0,   0,  0, 0,   0, 0,   0],\
                           [0, 0,   0,   0,  1, 0,   0, 0,   0],\
                           [0, 0,   0,   0,  0, 0,   0, 0,   1],\
                           [0, 0.5, 0, 0.5,  0, 0,   0, 0,   0],\
                           [0, 0, 0.5,   0,  0, 0, 0.5, 0,   0],\
                           [0, 0,   0,   0,  0, 0.5, 0, 0.5, 0] ], dtype = 'd')

        # Pseudo inverse from numpy: D*Dplus=Id6
	Dplus = numpy.matrix([ [ 1.,  0.,  0.,  0.,  0.,  0.],
		    	       [ 0.,  0.,  0.,  1.,  0.,  0.],
		    	       [ 0.,  0.,  0.,  0.,  1.,  0.],
		    	       [ 0.,  0.,  0.,  1.,  0.,  0.],
		    	       [ 0.,  1.,  0.,  0.,  0.,  0.],
		    	       [ 0.,  0.,  0.,  0.,  0.,  1.],
		    	       [ 0.,  0.,  0.,  0.,  1.,  0.],
		    	       [ 0.,  0.,  0.,  0.,  0.,  1.],
		    	       [ 0.,  0.,  1.,  0.,  0.,  0.]], dtype = 'd')
        #print 'D =\n%s' % D
        #print 'Dplus =\n%s' % Dplus
        #print 'D*Dplus =\n%s' % (D*Dplus)

	# Defining G+
        g1 = 0.5/(L2 - L1)
        g2 = 0.5/(L0 - L2)
        g3 = 0.5/(L1 - L0)
	Gplus = numpy.matrix([ [1, 0 , 0 , 0 , 0, 0 , 0 , 0 , 0],
                               [0, 0 , 0 , 0 , 1, 0 , 0 , 0 , 0],
                               [0, 0 , 0 , 0 , 0, 0 , 0 , 0 , 1],
                               [0, 0 , 0 , 0 , 0, g1, 0 , g1, 0],
                               [0, 0 , g2, 0 , 0, 0 , g2, 0 , 0],
                               [0, g3, 0 , g3, 0, 0 , 0 , 0 , 0] ],
                             dtype = 'd')
        #print 'Gplus =\n%s' % Gplus

	# Defining F^-1
        Finverse = Gplus * numpy.kron(S, S) * Dplus

	# Derive the error propagation matrix K
	K_low   = numpy.c_[M_IDENTITY_3_3, M_ZEROS_3_3]
        K_left = numpy.c_['0', numpy.c_['1', M_ZEROS_9_3, S_p], K_low]
        #print 'K_left =\n%s' % K_left
	K = K_left * Finverse
        #print 'K =\n%s' % K
        
	## Define the covariance matrix of the errors on the inertia tensor
        cIxx_xx = 0.0
	cIxx_yy = 0.0
	cIxx_zz = 0.0
	cIxx_xy = 0.0
	cIxx_xz = 0.0
	cIxx_yz = 0.0
	
        cIyy_yy = 0.0
        cIyy_zz = 0.0
        cIyy_xy = 0.0
        cIyy_xz = 0.0
        cIyy_yz = 0.0

	cIzz_zz = 0.0
	cIzz_xy = 0.0
	cIzz_xz = 0.0
	cIzz_yz = 0.0	

	cIxy_xy = 0.0
	cIxy_xz = 0.0
	cIxy_yz = 0.0

	cIxz_xz = 0.0
	cIxz_yz = 0.0

	cIyz_yz = 0.0
	
        for dataPoint in dataVec:
            hit = vector2point(dataPoint.getPoint()) - self.Centroid
            x = hit.x()
            y = hit.y()
            z = hit.z()
            w    = dataPoint.getWeight()
	    # Need something smarter, here!
	    dx = 5.0
	    dy = 5.0
	    dz = 5.0
	    dw = 0.1*w
	    
            # Ixx-Others
            cIxx_xx +=  4*w**2 * (y**2 * dy**2 + z**2 * dz**2) + (y**2 + z**2)**2  * dw**2 
            cIxx_yy +=  4*w**2 * z**2 * dz**2         + (y**2 + z**2)*(x**2 + z**2) * dw**2
            cIxx_zz +=  4*w**2 * y**2 * dy**2         + (y**2 + z**2)*(x**2 + y**2) * dw**2
            cIxx_xy += -2*w**2 * x*y  * dy**2         + (y**2 + z**2)*(-x*y)        * dw**2
            cIxx_xz += -2*w**2 * x*z  * dz**2         + (y**2 + z**2)*(-x*z)        * dw**2
            cIxx_yz += -2*w**2 * y*z  * (dy**2+dz**2) + (y**2 + z**2)*(-y*z)        * dw**2

            # Iyy-Others
            cIyy_yy +=  4*w**2 * (x**2 * dx**2 + z**2 * dz**2) + (x**2 + z**2)**2   * dw**2
            cIyy_zz +=  4*w**2 * x**2 * dx**2	      + (x**2 + z**2)*(x**2 + y**2) * dw**2    
            cIyy_xy += -2*w**2 * x*y  * dx**2	      + (x**2 + z**2)*(-x*y)	    * dw**2	   
            cIyy_xz += -2*w**2 * x*z  * (dx**2+dz**2) + (x**2 + z**2)*(-x*z)	    * dw**2	    
            cIyy_yz += -2*w**2 * y*z  * dz**2	      + (x**2 + z**2)*(-y*z)	    * dw**2

             # Iyy-Others
            cIzz_zz += 4*w**2 * (x**2 * dx**2 + y**2 * dy**2) + (x**2 + y**2)**2    * dw**2
	    cIzz_xy += -2*w**2 * x*y  * (dx**2+dy**2) + (x**2 + y**2)*(-x*y)	    * dw**2
            cIzz_xz += -2*w**2 * x*z  * dx**2	      + (x**2 + y**2)*(-x*z)	    * dw**2
            cIzz_yz += -2*w**2 * y*z  * dy**2	      + (x**2 + y**2)*(-y*z)	    * dw**2

            # Ixy-Others
            cIxy_xy += w**2 * (y**2 * dx**2 + x**2 * dy**2)   + x**2 * y**2         * dw**2
            cIxy_xz += w**2 * y*z * dx**2                     + x**2 * y*z          * dw**2
            cIxy_yz += w**2 * x*z * dy**2                     + x * y**2 * z        * dw**2

            # Ixz-Others
            cIxz_xz += w**2 * (z**2 * dx**2 + x**2 * dz**2)   + x**2 * z**2         * dw**2
	    cIxz_yz += w**2 * x*y * dz**2                     + x * y * z**2        * dw**2
	    
            # Iyz-Iyz
            cIyz_yz += w**2 * (y**2 * dz**2 + z**2 * dy**2)   + y**2 * z**2         * dw**2
	    

        self.VdICovMatrix = numpy.array([\
	     [cIxx_xx, cIxx_yy, cIxx_zz, cIxx_xy, cIxx_xz, cIxx_yz],
	     
             [cIxx_yy, cIyy_yy, cIyy_zz, cIyy_xy, cIyy_xz, cIyy_yz],
	     
             [cIxx_zz, cIyy_zz, cIzz_zz, cIzz_xy, cIzz_xz, cIzz_yz],
	     
	     [cIxx_xy, cIyy_xy, cIzz_xy, cIxy_xy, cIxy_xz, cIxy_yz],
	     
             [cIxx_xz, cIyy_xz, cIzz_xz, cIxy_xz, cIxz_xz, cIxz_yz],
	     
	     [cIxx_yz, cIyy_yz, cIzz_yz, cIxy_yz, cIxz_yz, cIyz_yz]],
             dtype = 'd')
        
	#print 'VdI covariance matrix'
	#print self.VdICovMatrix
	# Propagate errors, i.e.
	# Calculate the covariance matrix of Eigen values and vectors
	# Covariance matrix of:
        # [de0x, de0y, de0z, de1x, de1y, de1z, de2x, de2y, de2z, dL0, dL1, dL2]
	errorCovarianceMatrix = K * self.VdICovMatrix * K.transpose()
        #print 'Sigma(dvec(S), dLambda)\n%s' % errorCovarianceMatrix
        print 'dxdir = %s' % sqrt(errorCovarianceMatrix[1,1])
        print 'dydir = %s' % sqrt(errorCovarianceMatrix[4,4])
        print 'dzdir = %s' % sqrt(errorCovarianceMatrix[7,7])
        print 'dLambda1 = %.3e' % sqrt(errorCovarianceMatrix[9,9])
        print 'dLambda2 = %.3e' % sqrt(errorCovarianceMatrix[10,10])
        print 'dLambda3 = %.3e' % sqrt(errorCovarianceMatrix[11,11])
        return errorCovarianceMatrix

    def __str__(self):
        text = ''
        text += 'Centroid = %s\n' % self.Centroid
        text += 'Moments  = %s\n' % self.Moment
        for i in range(3):
            text += 'Axis %d   = %s\n' % (i, self.Axis[i])
        text += 'Longitudinal RMS = %.3e\n' % self.RmsLong
        text += 'Transverse RMS   = %.3e\n' % self.RmsTrans
        text += 'Long. asymmetry  = %f\n' % self.RmsLongAsym
        text += 'Long skewness    = %.3e' % self.SkewnessLong
        return text



class MomentsClusterInfo:

    def __init__(self, xtalCol, calRecon = None,
                 centrDistCut = 5.0*TOWER_PITCH, clipDataVec = False,
                 lastClipStep = False, pruneLong = False):
        self.RootCanvas = None
        self.XtalCol = xtalCol
        self.CalRecon = calRecon
        self.NotEnoughXtals = False
        iniCentroid = self.getInitialCentroid()
        dataVec = self.getCalMomentsDataVec(iniCentroid, centrDistCut)
        if len(dataVec) < MIN_NUM_XTALS:
            print 'Not enough xtals found (%d).' % len(dataVec)
            self.NotEnoughXtals = True
            self.MomentsAnalysis = CalMomentsAnalysis(self, clipDataVec)
            self.Centroid = self.MomentsAnalysis.estimateCentroid(dataVec)
            return
        self.MomentsAnalysis = CalMomentsAnalysis(self, clipDataVec)
        (chiSq, dataVec) =\
                self.MomentsAnalysis.doIterativeMomentsAnalysis(dataVec,
                                                                iniCentroid)
        self.Centroid = None
        self.Axis = None
        self.RmsLong = 0.
        self.RmsTrans = 0.
        self.RmsLongAsym = 0.
        self.SkewnessLong = -9999.
        if chiSq >= 0.:
            # Get info on iterations
            nIterations = self.MomentsAnalysis.getNumIterations()
            nDropped = self.MomentsAnalysis.getNumDroppedPoints()
            # Get the iterative moments centroid and axis from iterations
            self.Centroid = self.MomentsAnalysis.getMomentsCentroid()
            self.Axis = self.MomentsAnalysis.getMomentsAxis()
            self.Moment = copy.copy(self.MomentsAnalysis.Moment)
            # Recalculate the moments going back to using all the data points
            # but with the iterated moments centroid
            if nIterations > 1:
                #self.CovMatrix =\
                #     self.MomentsAnalysis.getCovarianceMatrix(dataVec)

                if lastClipStep:
                    # One more iteration with the clipped data vec.
                    clippedDataVec = self.getClippedDataVec(dataVec)
                    chiSq =\
                       self.MomentsAnalysis.doMomentsAnalysis(clippedDataVec,
                                                              self.Centroid)
                    self.MomentsAnalysis.NumIterations += 1
                    self.Centroid = self.MomentsAnalysis.getMomentsCentroid()
                    self.Axis = self.MomentsAnalysis.getMomentsAxis()
                    self.Moment = copy.copy(self.MomentsAnalysis.Moment)
                    # End of the new iteration.

                if pruneLong:
                    lastIter = self.MomentsAnalysis.IterationList[-1]
                    ppDataVec = lastIter.MomentsDataListPP
                    chiSq =\
                       self.MomentsAnalysis.doMomentsAnalysis(ppDataVec,
                                                              self.Centroid)
                    self.MomentsAnalysis.NumIterations += 1
                    self.Centroid = self.MomentsAnalysis.getMomentsCentroid()
                    self.Axis = self.MomentsAnalysis.getMomentsAxis()
                    self.Moment = copy.copy(self.MomentsAnalysis.Moment)

                dataVec = self.getCalMomentsDataVec(iniCentroid)
                chiSq = self.MomentsAnalysis.doMomentsAnalysis(dataVec,
                                                          self.Centroid)
            else:
                pass
                #self.CovMatrix =\
                #     self.MomentsAnalysis.getCovarianceMatrix(dataVec)
            # Extract the values for the moments with all hits present
            self.RmsLong = self.MomentsAnalysis.getLongitudinalRms()
            self.RmsTrans = self.MomentsAnalysis.getTransverseRms()
            self.RmsLongAsym = self.MomentsAnalysis.getLongAsymmetry()
            self.SkewnessLong = self.MomentsAnalysis.getLongSkewness()
        if self.CalRecon is not None:
            self.checkMomentsAnalysis()

    def cleanup(self):
        print 'Cleaning up...'
        for iter in self.MomentsAnalysis.IterationList:
            try:
                iter.Canvas.Close()
            except:
                pass
        try:
            self.RootCanvas.Close()
        except:
            pass

    def getInitialCentroid(self):
        # First estimation of the centroid.
        weightSum = 0.
        centroid = Point()
        for xtalData in self.XtalCol:
            dataPoint = CalMomentsData(xtalData)
            weight = dataPoint.getWeight()
            weightSum += weight
            centroid += vector2point(dataPoint.getPoint()) * weight
        centroid /= weightSum
        return centroid

    def getCalMomentsDataVec(self, centroid, maxCentrDist = 5*TOWER_PITCH,
                             minFracWeight = 0.0, preprocess = True):
        dataVec = []
        # This is a shortcut avoiding the initial pruning that Tracy does in
        # order to remove the outliers before the moments analysis.
        if not preprocess:
            for xtalData in self.XtalCol:
                dataVec.append(CalMomentsData(xtalData))
            return dataVec

        # Pre-calculate the weight sum (for the cut on minFracWeight).
        weightSum = 0.
        for xtalData in self.XtalCol:
            momentsData = CalMomentsData(xtalData)
            weightSum += momentsData.getWeight()
        minWeight = weightSum*minFracWeight
        
        # And this is the (almost) full code in the CalRecon.
        # Set the initial axis pointing up.
        axis = Point(0, 0, 1)
        # Go on with Tracy's code.
        rmsDist   = 0.
        weightSum = 0.
        for xtalData in self.XtalCol:
            momentsData = CalMomentsData(xtalData)
            # IMPORTANT!!
            # Portion of code handling the saturated crystals missing, here!
            distToAxis = momentsData.calcDistToAxis(centroid, axis)
            rmsDist   += momentsData.getWeight() * distToAxis * distToAxis
            centrDist = momentsData.getDistToPoint(centroid)
            weight = momentsData.getWeight()
            weightSum += weight
            # New if for the two constructor parameters.
            if (centrDist < maxCentrDist) and (weight > minWeight):
                dataVec.append(momentsData)
        if len(dataVec) < MIN_NUM_XTALS:
            return dataVec
        # Get the quick rms of crystals about fit axis
        rmsDist = sqrt(rmsDist / weightSum);
        # Use this to try to remove "isolated" crystals that might otherwise
        # pull the axis
        # Start by sorting the data vector according to distance to axis
        dataVec.sort()
        # Look for the point wherethere is a significant gap in distance to
        # the next xtal
        envelope = min(5.*rmsDist, TOWER_PITCH)
        lastDist = 0.
        checkit = len(dataVec)
        tempCounter = 0
        # Find the point in the vector where the distance to the axis
        # indicates an outlier
        for (i, momentsData) in enumerate(dataVec):
            dist = momentsData.getDistToAxis()
            distDiff = dist - envelope
            gapDist  = dist - lastDist
            if distDiff > 0. and gapDist > 3. * rmsDist:
                break
            lastDist = dist
            tempCounter += 1
        # Now remove the outlier crystals
        dataVec = dataVec[:tempCounter]
        return dataVec

    def getClippedDataVec(self, dataVec):
        clippedDataVec = []
        dx = self.MomentsAnalysis.Axis[1].x()
        dy = self.MomentsAnalysis.Axis[1].y()
        dz = self.MomentsAnalysis.Axis[1].z()
        xc = self.MomentsAnalysis.Centroid.x()
        yc = self.MomentsAnalysis.Centroid.y()
        zc = self.MomentsAnalysis.Centroid.z()
        xbot = xc - (zc - CAL_MODULE_CSI_BOTTOM_CORR)*(dx/dz)
        ybot = yc - (zc - CAL_MODULE_CSI_BOTTOM_CORR)*(dy/dz)
        zbot = CAL_MODULE_CSI_BOTTOM_CORR
        xtop = xc + (CAL_MODULE_CSI_TOP_CORR - zc)*(dx/dz)
        ytop = yc + (CAL_MODULE_CSI_TOP_CORR - zc)*(dy/dz)
        ztop = CAL_MODULE_CSI_TOP_CORR
        minCoord = -sqrt((xc - xbot)**2 + (yc - ybot)**2 + (zc - zbot)**2)
        maxCoord = sqrt((xc - xtop)**2 + (yc - ytop)**2 + (zc - ztop)**2)
        padding = 0.5*MOLIERE_RADIUS*sqrt(1/dz**2 - 1)
        minCoord += padding
        maxCoord -= padding
        for momentsData in dataVec:
            coord = momentsData.CoordAlongAxis
            if coord > minCoord and coord < maxCoord:
                clippedDataVec.append(momentsData)
        return clippedDataVec

    def compare(self, v1, v2, label, threshold = 1e-6):
        try:
            if abs(v1 - v2)/abs(v1 + v2) > threshold:
                print 'Difference in %s: [%s vs. %s]' % (label, v1, v2)
                return 1
            return 0
        except ZeroDivisionError:
            return 0

    def checkMomentsAnalysis(self):
        print 'Cross-checking moments analysis against the original values...'
        #cluster = self.CalRecon.getCalClusterCol().At(0)
        #cluster = self.XtalCol.At(0)
        cluster = self.CalRecon.getCalClusterCol().At(0)
        nDiff = 0
        nDiff += self.compare(self.RmsLong, cluster.getRmsLong(), 'RmsLong')
        nDiff += self.compare(self.RmsTrans, cluster.getRmsTrans(), 'RmsTrans')
        nDiff += self.compare(self.RmsLongAsym, cluster.getRmsLongAsym(),
                              'LongAsym')
        nDiff += self.compare(self.SkewnessLong, cluster.getSkewnessLong(),
                              'LongSkewness')
        reconCentroid = vector2point(cluster.getParams().getCentroid())
        nDiff += self.compare(self.Centroid.x(), reconCentroid.x(), 'Centr. x')
        nDiff += self.compare(self.Centroid.y(), reconCentroid.y(), 'Centr. y')
        nDiff += self.compare(self.Centroid.z(), reconCentroid.z(), 'Centr. z')
        reconAxis = vector2point(cluster.getParams().getAxis())
        nDiff += self.compare(self.Axis.x(), reconAxis.x(), 'Axis xdir')
        nDiff += self.compare(self.Axis.y(), reconAxis.y(), 'Axis ydir')
        nDiff += self.compare(self.Axis.z(), reconAxis.z(), 'Axis zdir')
        if nDiff:
            self.MomentsAnalysis.EqualToRecon = False
            print 'Done, %d difference(s) found.' % nDiff
            #print 'WARNING! Press enter to proceed, q to quit.'
            #if raw_input() == 'q':
            #    sys.exit()
        else:
            self.MomentsAnalysis.EqualToRecon = True
            print 'Done, no differences :-)'

    def drawLongProfile(self, expSkewness = None):
        if self.RootCanvas is None:
            self.RootCanvas = ROOT.TCanvas('cmoments', 'Moments analysis',
                                           1100, 600)
        self.RootCanvas.Clear()
        self.RootCanvas.Divide(3, 2)
        for (i, iteration) in enumerate(self.MomentsAnalysis.IterationList):
            self.RootCanvas.cd(i + 1)
            iteration.drawLongProfile(expSkewness)
        self.RootCanvas.cd()
        self.RootCanvas.Update()

    def __str__(self):
        text = ''
        text += 'Centroid         = %s\n' % self.Centroid
        text += 'Principal axis   = %s\n' % self.Axis
        text += 'Longitudinal RMS = %.3e\n' % self.RmsLong
        text += 'Transverse RMS   = %.3e\n' % self.RmsTrans
        text += 'Long. asymmetry  = %f\n' % self.RmsLongAsym
        text += 'Long skewness    = %.3e' % self.SkewnessLong
        return text
        


if __name__ == '__main__':
    from optparse import OptionParser
    usage = 'CalMomentsAnalysis.py reconFile <meritFile>'
    parser = OptionParser()
    parser.add_option('-c', '--skim-cut', type = str, dest = 'c',
                      default = '1',
                      help = 'a cut to filter the events')
    parser.add_option('-d', '--max-ctr-dist', type = float, dest = 'd',
                      default = 5.0,
                      help = 'pre-cut on the hit distance from the centroid')
    parser.add_option('-l', '--clip', action = 'store_true', dest = 'l',
                      default = False,
                      help = 'clip the data vec at each step')
    parser.add_option('-p', '--prune-long', action = 'store_true', dest = 'p',
                      default = False,
                      help = 'fix the xtals with bad long position')
    parser.add_option('-L', '--last-clip', action = 'store_true', dest = 'L',
                      default = False,
                      help = 'one more iteration after the last clip')
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
    reader = ReconReader(reconFilePath, meritFilePath, opts.c)
    eventNumber = 0
    answer = ''
    while answer != 'q':
        if reader.getEntry(eventNumber):
            calRecon = reader.getCalRecon()
            xtalCol = calRecon.getCalXtalRecCol()
            m = MomentsClusterInfo(xtalCol, calRecon, opts.d*TOWER_PITCH,
                                   opts.l, opts.L, opts.p)
            if not m.NotEnoughXtals:
                for iter in m.MomentsAnalysis.IterationList:
                    iter.draw(reader)
                    print '\n%s' % iter
                print 
                print '*** Final parameters after %d iteration(s):\n%s' %\
                      (m.MomentsAnalysis.getNumIterations(), m)
                m.drawLongProfile(reader.ExpectedSkewness)
            answer = raw_input('\nType q to quit, s to save, enter to go on')
            if answer == 's':
                for (i, iter) in enumerate(m.MomentsAnalysis.IterationList):
                    evtId = reader.getMeritVariable('EvtEventId')
                    filePath = 'Event_%s_iter%d' % (evtId, i)
                    iter.save(filePath)
            m.cleanup()
        eventNumber += 1
        

