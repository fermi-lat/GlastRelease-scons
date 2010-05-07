
import ROOT
import numpy
import random
import copy

from math import sqrt, sin, cos, atan, pi

RAD_TO_DEG = 180./pi
DEG_TO_RAD = pi/180.

ROOT.gStyle.SetCanvasColor(ROOT.kWhite)

X_LIST    = [5., 20., 34., 36., 44., 56., 78.]
Y_LIST    = [7., 10., 23., 20., 42., 40., 45.]
W_LIST    = [3., 10., 20., 8. , 8. , 23., 12.]
POS_ERR   = 1.0
W_REL_ERR = 0.0


D = numpy.array([ [1., 0. , 0. , 0.],
                  [0., 0. , 0. , 1.],
                  [0., 0.5, 0.5, 0.] ],
                dtype = 'd')
D_PLUS = numpy.array([ [1., 0., 0.],
                       [0., 0., 1.],
                       [0., 0., 1.],
                       [0., 1., 0.] ],
                     dtype = 'd')




class Point2D:
    
    def __init__(self, x, y):
        self.X = x
        self.Y = y

    def __add__(self, other):
        return Point2D(self.X + other.X, self.Y + other.Y)

    def __mul__(self, value):
        return Point2D(self.X*value, self.Y*value)

    def __div__(self, value):
        return Point2D(self.X/value, self.Y/value)
        
    def __str__(self):
        return '(%.3f, %.3f)' % (self.X, self.Y)



class Vector2D(Point2D):

    def __init__(self, xdir, ydir, origin):
        Point2D.__init__(self, xdir, ydir)
        self.Origin = origin

    def draw(self, scale = 1):
        self.Arrow = ROOT.TArrow(self.Origin.X, self.Origin.Y,
                                 self.Origin.X + scale*self.X,
                                 self.Origin.Y + scale*self.Y, 0.02, '|>')
        self.Arrow.Draw()



class Hit2D:

    def __init__(self, x, y, w, dx, dy, dw):
        self.X = x + random.gauss(0, dx)
        self.Y = y + random.gauss(0, dy)
        self.Weight = w + random.gauss(0, dw)
        self.ErrX = dx
        self.ErrY = dy
        self.ErrWeight = dw

    def draw(self):
        self.Marker = ROOT.TMarker(self.X, self.Y, 25)
        self.Marker.SetMarkerSize(0.1*self.Weight)
        self.Marker.Draw()

    def __str__(self):
        return 'x = %.3f +- %.3f, y = %.3f +- %.3f, w = %.3f +- %.3f' %\
            (self.X, self.ErrX, self.Y, self.ErrY, self.Weight, self.ErrWeight)


    
class Cluster2D:
    
    def __init__(self, hitList = None):
        if hitList is None:
            self.HitList = []
            for i in xrange(len(X_LIST)):
                hit = Hit2D(X_LIST[i], Y_LIST[i], W_LIST[i],
                            POS_ERR, POS_ERR, W_REL_ERR*W_LIST[i])
                self.HitList.append(hit)
        else:
            self.HitList = hitList
        self.__calcBaricenter()
        self.__calcInertiaTensor()
    
    def __calcBaricenter(self):
        self.Baricenter = Point2D(0, 0)
        self.WeightSum = 0
        for hit in self.HitList:
            x = hit.X
            y = hit.Y
            w = hit.Weight
            self.Baricenter += Point2D(x, y)*w
            self.WeightSum += w
        if self.WeightSum > 0:
            self.Baricenter /= self.WeightSum

    def __calcInertiaTensor(self):
        Ixx = 0.0
        Iyy = 0.0
        Ixy = 0.0
        cIxx_xx = 0.0
        cIyy_yy = 0.0
        cIxy_xy = 0.0
        cIxx_yy = 0.0
        cIxx_xy = 0.0
        cIxy_yy = 0.0
        for hit in self.HitList:
            x = hit.X - self.Baricenter.X
            y = hit.Y - self.Baricenter.Y
            w = hit.Weight
            dx = hit.ErrX
            dy = hit.ErrY
            dw = hit.ErrWeight
            Ixx += w * y**2
            Iyy += w * x**2
            Ixy -= w * x * y
            cIxx_xx += 4 * w**2 * y**2 * dy**2 + y**4 * dw**2
            cIyy_yy += 4 * w**2 * x**2 * dx**2 + x**4 * dw**2
            cIxy_xy += w**2 * y**2 * dx**2 + w**2 * x**2 * dy**2 +\
                x**2 * y**2 * dw**2
            cIxx_yy += x**2 * y**2 * dw**2
            cIxx_xy -= 2 * w**2 * x * y * dy**2 + x * y**3 * dw**2
            cIxy_yy -= 2 * w**2 * x * y * dx**2 + x**3 * y * dw**2
        self.I = numpy.array([ [Ixx, Ixy],
                               [Ixy, Iyy] ],
                             dtype = 'd')
        self.VecI = numpy.array([ [Ixx],
                                  [Ixy],
                                  [Ixy],
                                  [Iyy] ],
                                dtype = 'd')
        self.VdI = numpy.dot(D, self.VecI)
        self.VdICovMatrix = numpy.array([ [cIxx_xx, cIxx_yy, cIxx_xy],
                                          [cIxx_yy, cIyy_yy, cIxy_yy],
                                          [cIxx_xy, cIxy_yy, cIxy_xy] ],
                                        dtype = 'd')
        delta = sqrt((Ixx - Iyy)**2 + 4*Ixy**2)
        l0 = 0.5*(Ixx + Iyy - delta)
        l1 = 0.5*(Ixx + Iyy + delta)
        self.Lambda = numpy.array( [ [l0, 0.],
                                     [0., l1] ],
                                   dtype = 'd')
        self.VecLambda = numpy.array([ [l0],
                                       [0.],
                                       [0.],
                                       [l1] ],
                                     dtype = 'd')
        self.VdLambda = numpy.dot(D, self.VecLambda)
        self.Eigenvalues = numpy.array([ [l0],
                                         [l1] ],
                                       dtype = 'd')
        phi = -0.5*atan(2*Ixy/(Iyy - Ixx))
        e0x = cos(phi)
        e0y = sin(phi)
        e1x = -sin(phi)
        e1y = cos(phi)
        self.MajorAxis = Vector2D(e0x, e0y, self.Baricenter)
        self.MinorAxis = Vector2D(e1x, e1y, self.Baricenter)
        self.S = numpy.array([ [e0x, e0y],
                               [e1x, e1y] ],
                             dtype = 'd')
        self.T = numpy.kron(self.S, self.S)
        self.V = numpy.dot(numpy.dot(D, self.T), D_PLUS)
        self.G = numpy.array([ [1., 0., 0.],
                               [0., 0., (l1-l0)],
                               [0., 0., (l1-l0)],
                               [0., 1., 0.] ],
                             dtype = 'd')
        self.TT = numpy.kron(self.S.T, self.S.T)
        self.F = numpy.dot(numpy.dot(D, self.TT), self.G)
        self.GPlus = numpy.array([ [1., 0., 0., 0.],
                                   [0., 0., 0., 1.],
                                   [0., 0.5/(l1-l0), 0.5/(l1-l0), 0.] ],
                                 dtype = 'd')
        self.FInv = numpy.dot(numpy.dot(self.GPlus, self.T), D_PLUS)
        self.Phi = RAD_TO_DEG*phi
 
    def draw(self):
        self.Canvas = ROOT.TCanvas('c', 'c', 500, 500)
        self.Canvas.SetMargin(0.1, 0.1, 0.1, 0.1)
        self.BaseHist = ROOT.TH2F('hbase', 'cluster', 100, 0, 100, 100, 0, 100)
        self.BaseHist.Draw()
        for hit in self.HitList:
            hit.draw()
        self.BaricenterBox =\
            ROOT.TMarker(self.Baricenter.X, self.Baricenter.Y, 20)
        self.BaricenterBox.SetMarkerColor(ROOT.kRed)
        self.BaricenterBox.Draw()
        self.MajorAxis.draw(20)
        self.MinorAxis.draw(10)
        self.Canvas.Update()
    
    def __str__(self):
        text = 'Hits:\n'
        for hit in self.HitList:
            text += '%s\n' % hit
        text += 'Baricenter    : %s\n' % self.Baricenter
        text += 'Inertia tensor:\n%s\n' % self.I
        text += 'Inertia tensor covariance matrix:\n%s\n' % self.VdICovMatrix
        text += 'Eigenvalues   :\n%s\n' % self.Eigenvalues
        text += 'Rotation angle: %.3f deg\n' % self.Phi
        text += 'Major axis    : %s\n' % self.MajorAxis
        text += 'Minor axis    : %s\n' % self.MinorAxis
        return text



## Small test routines.

def mcVdICovMatrix(numTrials = 1000):
    meanIxx = 0.
    meanIyy = 0.
    meanIxy = 0.
    for i in xrange(numTrials):
        c = Cluster2D()
        meanIxx += c.I[0,0]
        meanIyy += c.I[1,1]
        meanIxy += c.I[0,1]
    meanIxx /= float(numTrials)
    meanIyy /= float(numTrials)
    meanIxy /= float(numTrials)
    dIxxdIxx = 0.
    dIxxdIyy = 0.
    dIxxdIxy = 0.
    dIyydIyy = 0.
    dIyydIxy = 0.
    dIxydIxy = 0.
    for i in xrange(numTrials):
        c = Cluster2D()
        Ixx = c.I[0,0]
        Iyy = c.I[1,1]
        Ixy = c.I[0,1]
        dIxxdIxx += (Ixx - meanIxx)*(Ixx - meanIxx)
        dIxxdIyy += (Ixx - meanIxx)*(Iyy - meanIyy) 
        dIxxdIxy += (Ixx - meanIxx)*(Ixy - meanIxy) 
        dIyydIyy += (Iyy - meanIyy)*(Iyy - meanIyy)
        dIyydIxy += (Iyy - meanIyy)*(Ixy - meanIxy)
        dIxydIxy += (Ixy - meanIxy)*(Ixy - meanIxy)
    dIxxdIxx /= float(numTrials)
    dIxxdIyy /= float(numTrials)
    dIxxdIxy /= float(numTrials)
    dIyydIyy /= float(numTrials)
    dIyydIxy /= float(numTrials)
    dIxydIxy /= float(numTrials)
    m = numpy.array([ [dIxxdIxx, dIxxdIyy, dIxxdIxy],
                      [dIxxdIyy, dIyydIyy, dIyydIxy],
                      [dIxxdIxy, dIyydIxy, dIxydIxy] ],
                    dtype = 'd')
    return m


def testPaperEqs(cluster):
    print '*** Testing D...'
    print 'D+ D =\n%s' % numpy.dot(D_PLUS, D)
    print
    print '*** Testing equation (1), (18), (19)...'
    print 'I =\n%s' % c.I
    print 'Lambda =\n%s' % c.Lambda
    print 'S I S^T =\n%s' %  (numpy.dot(c.S, numpy.dot(c.I, c.S.T)))
    print 'S^T Lambda S =\n%s' % (numpy.dot(c.S.T, numpy.dot(c.Lambda, c.S)))
    print
    print '*** Testing equation (3)...'
    print 'vec(I) =\n%s' % c.VecI
    print
    print '*** Testing equation (4)...'
    print 'vec(Lambda) =\n%s' % c.VecLambda
    print 'T vec(I) =\n%s' % (numpy.dot(c.T, c.VecI))
    print
    print '*** Testing equation (10)...'
    print 'v_d(I) =\n%s' % c.VdI
    print 'D vec(I) =\n%s' % numpy.dot(D, c.VecI)
    print 
    print '*** Testing equation (12)...'
    print 'v_d(Lambda) =\n%s' % c.VdLambda
    print 'V v_d(I) =\n%s' % numpy.dot(c.V, c.VdI)
    print
    print '*** Testing equation (15)...'
    print 'Sigma(v_d(I)) =\n%s' % c.VdICovMatrix
    s = numpy.dot(numpy.dot(c.V, c.VdICovMatrix), c.V.T)
    print 'V v_d(Lambda) V^T =\n%s' % s
    print 'Diagonal terms (square roots): (%.3f, %.3f, %.3f)' %\
          (sqrt(s[0,0]), sqrt(s[1,1]), sqrt(s[2,2]))
    print 
    print 'Testing equation 41...'
    print 'G+ G =\n%s' % numpy.dot(c.GPlus, c.G)
    print
    print 'Testing F inverse...'
    print 'F-1 F =\n%s' % numpy.dot(c.FInv, c.F)
    print
    print 'Testing equation (45)...'
    s = numpy.dot(numpy.dot(c.FInv, c.VdICovMatrix), c.FInv.T)
    print 'Finv VdICovMatrix Finv^T =\n%s' % s
    print 'Diagonal terms (square roots): (%.3f, %.3f, %.3f)' %\
          (sqrt(s[0,0]), sqrt(s[1,1]), RAD_TO_DEG*sqrt(s[2,2]))


if __name__ == '__main__':
    c = Cluster2D()
    print c
    c.draw()
    testPaperEqs(c)

    # The following done with 100,000 events gives:
    #[[  654306.79366708     2046.95299436  -444317.33476658]
    # [    2046.95299436  1410589.74012409  -445804.20738293]
    # [ -444317.33476658  -445804.20738293   514830.77784087]]
    #
    #m = mcVdICovMatrix(100000)
    #print 'Mc I covariance matrix:\n%s\n' % m

    
