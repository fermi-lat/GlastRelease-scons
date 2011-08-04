from math import *

import random
import numpy

import ROOT
ROOT.gROOT.SetStyle('Plain')



D_THETA = radians(2.0)
D_PHI   = radians(2.0)
Z_INT   = 1.0
NUM_EVT = 100000

# Basic formulas:
#
# xdir = sin(theta)cos(phi)
# ydir = sin(theta)sin(phi)
# zdir = cos(theta)
#


def getIntPoint(theta, phi):
    """ Return the intersection point (x, y) with the plane at z = 1 for a
    given direction.
    """
    x = Z_INT*tan(theta)*cos(phi)
    y = Z_INT*tan(theta)*sin(phi)
    return x, y

def getCovMatrix(theta, phi, dtheta, dphi):
    """ Return the covariance matrix for a given set of directions/errors.
    """
    st = sin(theta)
    ct = cos(theta)
    sp = sin(phi)
    cp = cos(phi)
    st2 = st**2.
    ct2 = ct**2.
    sp2 = sp**2.
    cp2 = cp**2.
    dt2 = dtheta**2.
    dp2 = dphi**2.
    m = numpy.matrix([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
    m[0, 0] = ct2*cp2*dt2 + st2*sp2*dp2
    m[1, 1] = ct2*sp2*dt2 + st2*cp2*dp2
    m[2, 2] = st2*dt2 
    m[0, 1] = ct2*cp*sp*dt2 - st2*cp*sp*dp2
    m[0, 2] = -st*ct*cp*dt2
    m[1, 2] = -st*ct*sp*dt2
    m[1, 0] = m[0, 1]
    m[2, 0] = m[0, 2]
    m[2, 1] = m[1, 2]
    return m


hresx = ROOT.TH1F('hresx', 'hresx', 100, -5, 5)
hresy = ROOT.TH1F('hresy', 'hresy', 100, -5, 5)

for i in xrange(NUM_EVT):
    u = random.uniform(0, 1)
    v = random.uniform(0, 1)
    theta0 = 2*pi*u
    phi0 = acos(2*v - 1)
    (x0, y0) = getIntPoint(theta0, phi0)
    theta = random.gauss(theta0, D_THETA)
    phi = random.gauss(phi0, D_PHI)
    (x, y) = getIntPoint(theta, phi)
    xdir = sin(theta)*cos(phi)
    ydir = sin(theta)*sin(phi)
    zdir = cos(theta)
    m = getCovMatrix(theta, phi, D_THETA, D_PHI)

    Dx = Z_INT/zdir
    Dz = -xdir*Z_INT/zdir**2.
    dx2  = Dx*Dx*m[0, 0] + Dz*Dz*m[2, 2] + 2*Dx*Dz*m[0, 2]
    dx   = sqrt(dx2)
    resx = (x - x0)/dx
    hresx.Fill(resx)

    Dy = Z_INT/zdir
    Dz = -ydir*Z_INT/zdir**2.
    dy2  = Dy*Dy*m[1, 1] + Dz*Dz*m[2, 2] + 2*Dy*Dz*m[1, 2]
    dy   = sqrt(dy2)
    resy = (y - y0)/dy
    hresy.Fill(resy)


c = ROOT.TCanvas('c', 'c', 1200, 600)
c.Divide(2, 1)
c.cd(1); hresx.Draw(); hresx.Fit('gaus')
c.cd(2); hresy.Draw(); hresy.Fit('gaus')
c.Update()
