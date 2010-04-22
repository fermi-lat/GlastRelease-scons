
# Author:
# Luca Baldini (luca.baldini@pi.infn.it)
#
# Purpose:
# Test program for the skewness in the moments analysis
# (see CalRecon/doc/moments_analysis.pdf)


# From Review of Particle Physics, July 2004, pag 250.
#
# dE/dt = E_0 * b * (b*t)^a*exp(-b*t) / (Gamma(a+1))
# par     0     1    1 x  2      1             2
# t = x/X0
#
# t_max = (a-1)/b = ln(y) + C
# y = E/Ec
# C = -0.5 for electrons
# C = +0.5 for photons 
#
# Following values are for CsI (from Review of Particle Physics, July 2004,
# pag 98).
# Radiation length = 8.39  g/cm^2 = 1.86 cm
# Critical energy  = 11.17 MeV for electrons (10.80 MeV for positrons)
# Moliere radius   = 15.92 g/cm^2 = 3.53 cm


from optparse import OptionParser

parser = OptionParser()
parser.add_option('-e', '--energy', type = float, dest = 'e', default = 1,
                  help = 'the initial particle energy (in GeV)')
parser.add_option('-m', '--tmin', type = float, dest = 'm', default = 0.0,
                  help = 'the minimum t for sampling (in X0)')
parser.add_option('-M', '--tmax', type = float, dest = 'M', default = 50.0,
                  help = 'the maximum t for sampling (in X0)')
parser.add_option('-p', '--particle-type', type = str, dest = 'p',
                  default = 'e',
                  help = 'the particle type (e or gamma)')
(opts, args) = parser.parse_args()

import ROOT
import sys

from math import sqrt, log, exp


INITIAL_ENERGY = float(opts.e)
T_MIN          = float(opts.m)
T_MAX          = float(opts.M)
if opts.p == 'e':
    C = -0.5
elif opts.p == 'gamma':
    C = 0.5
else:
    sys.exit('Unkown particle type (%s). Abort.' % opts.p)
B = 0.5

print 'Initializing for %s at %.3f GeV (tmin = %.3f, tmax = %.3f)' %\
      (opts.p, INITIAL_ENERGY, T_MIN, T_MAX)
    
CRITICAL_ENERGY = 0.01117 # GeV
GAMMA_FORMULA   = '[0]*x**[2]*exp(-[1]*x)'


y     = INITIAL_ENERGY/CRITICAL_ENERGY
tpeak = log(y) + C
b     = B
a     = b*tpeak


p  = ROOT.TF1('p' , GAMMA_FORMULA, T_MIN, T_MAX)
p1 = ROOT.TF1('p1', '(x)*%s' % GAMMA_FORMULA, T_MIN, T_MAX)
p2 = ROOT.TF1('p2', '(x**2)*%s' % GAMMA_FORMULA, T_MIN, T_MAX)
p3 = ROOT.TF1('p3', '(x**3)*%s' % GAMMA_FORMULA, T_MIN, T_MAX)
for function in [p, p1, p2, p3]:
    function.SetNpx(1000)
    function.SetParameter(0, 1)
    function.SetParameter(1, b)
    function.SetParameter(2, a)

norm   = p.Integral(T_MIN, T_MAX)
kf     = 1./norm
expt1  = p1.Integral(T_MIN, T_MAX)/norm
expt2  = p2.Integral(T_MIN, T_MAX)/norm
expt3  = p3.Integral(T_MIN, T_MAX)/norm
mu     = expt1
sigma2 = expt2 - expt1**2
sigma  = sqrt(sigma2)
gamma  = (expt3 - 3*mu*sigma2 - mu**3)/(sigma**3)

print '*** t_1 = %.3f, t_2 = %.3f' % (T_MIN, T_MAX)
print '*** a = %.3f, b = %.3f, shower max. at %.3f X0' % (a, b, tpeak)
print '*** k_f = %.3f' % kf
print '*** mu = %.3f' % mu
print '*** sigma = %.3f' % sigma
print '*** gamma = %.3f' % gamma
print


# Integrate to get the expectation value of t^n with arbitrary center and
# order.
#
# Uses the trapezium rule with n steps.
def integrateMoment(order = 3, center = 0, nsteps = 10):
    tmin   = T_MIN
    tmax   = T_MAX
    energy = INITIAL_ENERGY
    b      = B
    c      = C 
    a      = b*(log(energy/CRITICAL_ENERGY) + c)
    t1     = tmin
    numSum = 0.0
    denSum = 0.0
    for i in xrange(nsteps):
        t2 = tmin + (i + 1)*(tmax - tmin)/nsteps
        p1 = (t1**a)*exp(-b*t1)
        p2 = (t2**a)*exp(-b*t2)
        numSum += 0.5*( p1*((t1 - center)**order) +\
                        p2*((t2 - center)**order) )*(t2 - t1)
        denSum += 0.5*(p1 + p2)*(t2 - t1)
        t1 = t2
    return numSum/denSum


# Numerical integration to get the skewness.
#
# Uses the trapezium rule with n steps, meant to test the code to be put
# in the AnalysisTuple.
def integrateSkewness(nsteps = 10):
    tmin   = T_MIN
    tmax   = T_MAX
    energy = INITIAL_ENERGY
    b      = B
    c      = C
    a      = b*(log(energy/CRITICAL_ENERGY) + c)
    t1     = tmin
    norm   = 0.0
    mom1   = 0.0
    mom2   = 0.0
    mom3   = 0.0
    stepSize = (tmax - tmin)/nsteps
    for i in xrange(nsteps):
        t2 = tmin + (i + 1)*stepSize
        p1 = (t1**a)*exp(-b*t1)
        p2 = (t2**a)*exp(-b*t2)
        norm += 0.5*(p1 + p2)
        mom1 += 0.5*(p1*t1 + p2*t2)
        mom2 += 0.5*(p1*t1*t1 + p2*t2*t2)
        mom3 += 0.5*(p1*t1*t1*t1 + p2*t2*t2*t2)
        t1 = t2
    norm *= stepSize
    mom1 *= stepSize
    mom2 *= stepSize
    mom3 *= stepSize
    mom1 /= norm
    mom2 /= norm
    mom3 /= norm
    mom2 = mom2 - mom1*mom1
    gamma = (mom3 - 3*mom1*mom2 - mom1*mom1*mom1)/(mom2**1.5)
    return gamma


# Numerical integration to get the skewness.
#
# Uses a different method (discrete sampling), meant to test the code to be
# put in the AnalysisTuple.
def integrateSkewness2(nsteps = 8):
    tmin   = T_MIN
    tmax   = T_MAX
    energy = INITIAL_ENERGY
    b      = B
    c      = C
    a      = b*(log(energy/CRITICAL_ENERGY) + c)
    t1     = tmin
    norm   = 0.0
    mom1   = 0.0
    mom2   = 0.0
    mom3   = 0.0
    stepSize = (tmax - tmin)/nsteps
    for i in xrange(nsteps):
        t2 = tmin + (i + 1)*stepSize
        t = 0.5*(t1 + t2)
        p = (t**a)*exp(-b*t)
        norm += p
        mom1 += p*t
        mom2 += p*t*t
        mom3 += p*t*t*t
        t1 = t2
    norm *= stepSize
    mom1 *= stepSize
    mom2 *= stepSize
    mom3 *= stepSize
    mom1 /= norm
    mom2 /= norm
    mom3 /= norm
    mom2 = mom2 - mom1*mom1
    gamma = (mom3 - 3*mom1*mom2 - mom1*mom1*mom1)/(mom2**1.5)
    return gamma


# Numerical integration to get the skewness.
#
# Uses the trapezium rule with n steps, meant to test the code to be put
# in the AnalysisTuple.
def integrateSkewness3(nsteps = 10):
    tmin   = T_MIN
    tmax   = T_MAX
    energy = INITIAL_ENERGY
    b      = B
    c      = C
    a      = b*(log(energy/CRITICAL_ENERGY) + c)
    t1     = tmin
    norm   = 0.0
    mom1   = 0.0
    mom2   = 0.0
    mom3   = 0.0
    stepSize = (tmax - tmin)/nsteps
    for i in xrange(nsteps):
        t2 = tmin + (i + 1)*stepSize
        p1 = (t1**a)*exp(-b*t1)
        p2 = (t2**a)*exp(-b*t2)
        t  = 0.5*(t1 + t2)
        p  = 0.5*(p1 + p2)
        norm += p
        mom1 += p*t
        mom2 += p*t*t
        mom3 += p*t*t*t
        t1 = t2
    norm *= stepSize
    mom1 *= stepSize
    mom2 *= stepSize
    mom3 *= stepSize
    mom1 /= norm
    mom2 /= norm
    mom3 /= norm
    mom2 = mom2 - mom1*mom1
    gamma = (mom3 - 3*mom1*mom2 - mom1*mom1*mom1)/(mom2**1.5)
    return gamma


# Test the equations in the memo
# (see CalRecon/doc/moments_analysis.pdf)
def teq(eq):
    print '** Testing eq. || %s ||' % eq
    (left, right) = eq.split('=')
    left = left.strip()
    right = right.strip()
    leftv = eval(left)
    rightv = eval(right)
    print '%s = %.3f' % (left, leftv)
    print '%s = %.3f' % (right, rightv)
    print 'Fractional difference = %.2f' % (100*abs(leftv/rightv - 1))
    if abs(leftv/rightv - 1) < 0.001:
        print '----> OK'
    else:
        print '----> ERROR!'
    print


if T_MIN == 0 and T_MAX > 50:
    # Test exect formula (50 is essentially like infinity).
    teq('mu     = (a + 1)/b')
    teq('expt2  = (a + 2)*(a + 1)/b**2')
    teq('expt3  = (a + 3)*(a + 2)*(a + 1)/b**3')
    teq('sigma2 = (a + 1)/b**2')
    teq('gamma  = 2./sqrt(a + 1)')
    teq('expt2  = (a + 2)/b*expt1')
    teq('sigma2 = (a + 2)/b*mu - mu**2')
    teq('expt3  = (a + 3)/b*expt2')
    teq('gamma  = mu/(sigma**3)*((a + 3)*(a + 2)/b**2 - 3*sigma2 - mu**2)')
else:
    # Test integration by parts without dropping the first term.
    c1 = T_MAX**(a+2)*exp(-b*T_MAX)
    c2 = T_MIN**(a+2)*exp(-b*T_MIN)
    c  = (c1 - c2)/(b*norm)
    teq('expt2  = (a + 2)/b*expt1 - %s' % c)
    c1 = T_MAX**(a+3)*exp(-b*T_MAX)
    c2 = T_MIN**(a+3)*exp(-b*T_MIN)
    c  = (c1 - c2)/(b*norm)
    teq('expt3  = (a + 3)/b*expt2 - %s' % c)

m1Num = integrateMoment(1)
print 'Numerical mean: %.3f' % m1Num
m2Num = integrateMoment(2, m1Num)
print 'Numerical sigma**2: %.3f' % m2Num
gammaNum = integrateMoment(3, m1Num)/(m2Num**1.5)
for numSteps in [100, 50, 25, 10, 8, 7, 6, 5]:
    print 'Numerical skewness (%d steps): %.3f' %\
          (numSteps, integrateSkewness(numSteps))
    print 'Numerical skewness 2 (%d steps): %.3f' %\
          (numSteps, integrateSkewness2(numSteps))
    print 'Numerical skewness 3 (%d steps): %.3f' %\
          (numSteps, integrateSkewness3(numSteps))
print 


#p.Draw()
