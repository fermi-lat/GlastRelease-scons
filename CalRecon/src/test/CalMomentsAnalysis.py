
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

from CalLayout import *


ROOT.gStyle.SetCanvasColor(ROOT.kWhite)

M_PI = pi
CAL_TOWER_PITCH = 374.5

LIBRARIES = ['libcommonRootData.so', 'libreconRootData.so']
MERIT_VARS = ['EvtEventId',
              'McEnergy', 'CalEnergyRaw', 'CTBBestEnergy', 'CalEnergyCorr',
              'McXDir', 'McYDir', 'McZDir',
              'Tkr1XDir', 'Tkr1YDir', 'Tkr1ZDir',
              'Tkr1X0', 'Tkr1Y0', 'Tkr1Z0',
              'CalCsIRLn', 'CalLATRLn', 'TkrNumTracks']
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
        return info

    def getEntry(self, i):
        print 'ReconReader retriving event %d...' % i
        self.ReconChain.GetEvent(i)
        if self.MeritChain is not None:
            self.MeritChain.GetEntry(i)
            energy = self.getMeritVariable('CTBBestEnergy')
            tmax = self.getMeritVariable('CalLATRLn')
            tmin = tmax - self.getMeritVariable('CalCsIRLn')
            print energy, tmin, tmax
            (mean, sigma2, skewness) =\
                   getExpectedLongParameters(energy, tmin, tmax)
            self.ExpectedSkewness = skewness
            return self.MeritTreeFormula.EvalInstance()
        return 1

    def getCalRecon(self):
        return self.ReconEvent.getCalRecon()

    def getCalXtalRecCol(self):
        return self.getCalRecon().getCalXtalRecCol()



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
        self.Point = vector2point(self.getPoint())
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
            self.XZMarker.SetMarkerStyle(7)
            self.YZMarker.SetMarkerStyle(7)

    def getPoint(self):
        return self.XtalData.getPosition()

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
        for (i, dataPoint) in enumerate(dataVec):
            self.LongProfile.SetPoint(i, dataPoint.CoordAlongAxis,
                                      dataPoint.getWeight())
            dataPoint.setMaxWeight(momentsAnalysis.MaxWeight)
            self.MomentsDataList.append(copy.deepcopy(dataPoint))

    def draw(self, reconReader):
        dx = self.Axis[1].x()
        dy = self.Axis[1].y()
        dz = self.Axis[1].z()
        xc = self.Centroid.x()
        yc = self.Centroid.y()
        zc = self.Centroid.z()
        if reconReader.MeritChain is not None:
            tdx = reconReader.getMeritVariable('Tkr1XDir')
            tdy = reconReader.getMeritVariable('Tkr1YDir')
            tdz = reconReader.getMeritVariable('Tkr1ZDir')
            txc = reconReader.getMeritVariable('Tkr1X0')
            tyc = reconReader.getMeritVariable('Tkr1Y0')
            tzc = reconReader.getMeritVariable('Tkr1Z0')
        cName = 'cMomIter%d' % self.IterationNumber
        cTitle = 'Moments analysis---iteration %d' % self.IterationNumber
        self.Canvas = getCanvas(cName, cTitle)
        self.Canvas.cd(1)
        # Draw XZ view.
        CAL_LAYOUT.draw('xz')
        for dataPoint in  self.MomentsDataList:
            dataPoint.XZMarker.Draw()
        self.XZCentrMarker = ROOT.TMarker(xc, zc, 20)
        self.XZCentrMarker.SetMarkerColor(ROOT.kRed)
        self.XZCentrMarker.Draw()
        x1 = xc - (zc - Y_MIN)*(dx/dz)
        z1 = Y_MIN
        x2 = xc + (Y_MAX - zc)*(dx/dz)
        z2 = Y_MAX
        self.XZDirection = ROOT.TLine(x1, z1, x2, z2)
        self.XZDirection.SetLineColor(ROOT.kRed)
        self.XZDirection.SetLineWidth(1)
        self.XZDirection.Draw()
        if reconReader.MeritChain is not None:
            tx1 = txc - (tzc - Y_MIN)*(tdx/tdz)
            tz1 = Y_MIN
            tx2 = txc + (Y_MAX - tzc)*(tdx/tdz)
            tz2 = Y_MAX
            self.XZTkrDirection = ROOT.TLine(tx1, tz1, tx2, tz2)
            self.XZTkrDirection.SetLineColor(ROOT.kBlue)
            self.XZTkrDirection.SetLineWidth(1)
            self.XZTkrDirection.SetLineStyle(7)
            self.XZTkrDirection.Draw()
        self.Canvas.cd(2)
        # Draw YZ view.
        CAL_LAYOUT.draw('yz')
        for dataPoint in  self.MomentsDataList:
            dataPoint.YZMarker.Draw()
        self.YZCentrMarker = ROOT.TMarker(yc, zc, 20)
        self.YZCentrMarker.SetMarkerColor(ROOT.kRed)
        self.YZCentrMarker.Draw()
        y1 = yc - (zc - Y_MIN)*(dy/dz)
        y2 = yc + (Y_MAX - zc)*(dy/dz)
        self.YZDirection = ROOT.TLine(y1, z1, y2, z2)
        self.YZDirection.SetLineColor(ROOT.kRed)
        self.YZDirection.SetLineWidth(1)
        self.YZDirection.Draw()
        if reconReader.MeritChain is not None:
            ty1 = tyc - (tzc - Y_MIN)*(tdy/tdz)
            ty2 = tyc + (Y_MAX - tzc)*(tdy/tdz)
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
        text += 'Long skewness    = %.3e' % self.SkewnessLong
        return text

    

class CalMomentsAnalysis:

    def __init__(self):
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
        return chiSq

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

    def __init__(self, calRecon):
        self.RootCanvas = None
        self.CalRecon = calRecon
        self.NotEnoughXtals = False
        iniCentroid = self.getInitialCentroid()
        dataVec = self.getCalMomentsDataVec(iniCentroid)
        if len(dataVec) < MIN_NUM_XTALS:
            print 'Not enough xtals found (%d).' % len(dataVec)
            self.NotEnoughXtals = True
            return
        self.MomentsAnalysis = CalMomentsAnalysis()
        chiSq = self.MomentsAnalysis.doIterativeMomentsAnalysis(dataVec,
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
                self.CovMatrix = self.MomentsAnalysis.getCovarianceMatrix(dataVec)
                dataVec = self.getCalMomentsDataVec(iniCentroid)
                chiSq = self.MomentsAnalysis.doMomentsAnalysis(dataVec,
                                                          self.Centroid)
            else:
                pass
                self.CovMatrix = self.MomentsAnalysis.getCovarianceMatrix(dataVec)
            # Extract the values for the moments with all hits present
            self.RmsLong = self.MomentsAnalysis.getLongitudinalRms()
            self.RmsTrans = self.MomentsAnalysis.getTransverseRms()
            self.RmsLongAsym = self.MomentsAnalysis.getLongAsymmetry()
            self.SkewnessLong = self.MomentsAnalysis.getLongSkewness()
        self.checkMomentsAnalysis()
        

    def getInitialCentroid(self):
        # First estimation of the centroid.
        weightSum = 0.
        centroid = Point()
        for xtalData in self.CalRecon.getCalXtalRecCol():
            dataPoint = CalMomentsData(xtalData)
            weight = dataPoint.getWeight()
            weightSum += weight
            centroid += vector2point(dataPoint.getPoint()) * weight
        centroid /= weightSum
        return centroid

    def getCalMomentsDataVec(self, centroid, preprocess = True):
        dataVec = []
        # This is a shortcut avoiding the initial pruning that Tracy does in
        # order to remove the outlers before the moments analysis.
        if not preprocess:
            for xtalData in self.CalRecon.getCalXtalRecCol():
                dataVec.append(CalMomentsData(xtalData))
            return dataVec
        # And this is the (almost) full code in the CalRecon.
        # Set the initial axis pointing up.
        axis = Point(0, 0, 1)
        # Go on with Tracy's code.
        rmsDist   = 0.
        weightSum = 0.
        for xtalData in self.CalRecon.getCalXtalRecCol():
            momentsData = CalMomentsData(xtalData)
            # IMPORTANT!!
            # Portion of code handling the saturated crystals missing, here!
            distToAxis = momentsData.calcDistToAxis(centroid, axis)
            rmsDist   += momentsData.getWeight() * distToAxis * distToAxis
            weightSum += momentsData.getWeight()
            dataVec.append(momentsData)
        if len(dataVec) < MIN_NUM_XTALS:
            return dataVec
        # Get the quick rms of crystals about fit axis
        rmsDist = sqrt(rmsDist / weightSum);
        # Use this to try to remove "isolated" crystals that might otherwise
        # pull the axis
        # Start by sorting the data vector according to distance to axis
        dataVec.sort()
        # Look for the point where there is a significant gap in distance to
        # the next xtal
        envelope = min(5.*rmsDist, CAL_TOWER_PITCH)
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
        # Make calculated quantities class members.
        
        return dataVec

    def getFirstCluster(self):
        return self.CalRecon.getCalClusterCol()[0]

    def getUberCluster(self):
        # To be implemented.
        pass

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
            print 'Done, %d difference(s) found.' % nDiff
            print 'WARNING! Press enter to proceed.'
            raw_input()
        else:
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
                  default = '1', help = 'a cut to filter the events')
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
            m = MomentsClusterInfo(calRecon)
            if not m.NotEnoughXtals:
                for iter in m.MomentsAnalysis.IterationList:
                    print '\n%s' % iter
                    iter.draw(reader)
                print 
                print '*** Final parameters after %d iteration(s):\n%s' %\
                      (m.MomentsAnalysis.getNumIterations(), m)
                print '\n%s' % reader.getEventInfo()
                m.drawLongProfile(reader.ExpectedSkewness)
            answer = raw_input('\nType q to quit, enter to continue.')
        eventNumber += 1
        

