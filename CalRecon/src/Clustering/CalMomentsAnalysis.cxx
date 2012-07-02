/**
 *
 * File and Version Information:
 *      $Header$
 *
 */

#include "src/Clustering/CalMomentsAnalysis.h"

#include <algorithm>


void CalMomentsAnalysis::clear()
{
  // TBD: Need to think about what are the most sensible default values here!
  m_weightSum        = 0.;
  m_centroid         = Point(0.,0.,0.); 
  m_centroidErr      = CLHEP::HepMatrix(3,3,1);
  m_moment           = Vector(0.,0.,0.);
  m_axis[0]          = Vector(0.,0.,0.);
  m_axis[1]          = Vector(0.,0.,1.);
  m_axis[2]          = Vector(0.,0.,0.);
  m_axisErr          = CLHEP::HepMatrix(3,3,1);
  m_fullLength       = 0.;
  m_longRms          = 0.;
  m_transRms         = 0.;
  m_longRmsAsym      = 0.;
  m_longSkewness     = 0.;
  m_coreEnergyFrac   = 0;
  m_numIterations    = 0;
  m_numDroppedPoints = 0;
}

double CalMomentsAnalysis::doMomentsAnalysis(CalMomentsDataVec& dataVec,
                                             const Point& iniCentroid,
                                             double coreRadius)
{
  double chisq = -1.;
  
  // Check that we have enough points to proceed - need at least two.
  if ( dataVec.size() < 2 ) return chisq;
  
  // Initialize some local variables.
  double Ixx = 0.;
  double Iyy = 0.;
  double Izz = 0.;
  double Ixy = 0.;
  double Ixz = 0.;
  double Iyz = 0.;
  double weightSum = 0.;
  Point  centroid(0.,0.,0.);

  // Convert the core radius (passed as a parameter) in physical units (mm).
  // The Moliere radius for the CsI(Tl) is taken from the PDG Review of Particle Physics
  // (2008), table 28.4, pag. 288.
  coreRadius *= 35.7;
 
  // Loop through the data points a first time in order to determine
  // energy, centroid, and inertia tensor.
  CalMomentsDataVec::iterator vecIter;
  for(vecIter = dataVec.begin(); vecIter != dataVec.end(); vecIter++)
    {
      // Construct elements of (symmetric) "Inertia" Tensor:
      // See Goldstein, 1965, Chapter 5 (especially, eqs. 5-6, 7, 22, 26).
      // Analysis easy when translated to energy centroid.
 
      // Get pointer to the reconstructed data for given crystal.
      const CalMomentsData& dataPoint = *vecIter;
      double weight = dataPoint.getWeight();
      Vector hit = dataPoint.getPosition() - iniCentroid;
      double Rsq  = hit.mag2();
      double xprm = hit.x();
      double yprm = hit.y();
      double zprm = hit.z();

      // Update elements of the inertia tensor...
      Ixx += (Rsq - xprm*xprm) * weight;
      Iyy += (Rsq - yprm*yprm) * weight;
      Izz += (Rsq - zprm*zprm) * weight;
      Ixy -= xprm*yprm * weight;
      Ixz -= xprm*zprm * weight;
      Iyz -= yprm*zprm * weight;
      
      // ...and centroid/energy.
      weightSum += weight;
      centroid  += weight * dataPoint.getPosition();
    }

  // Render determinant of Inertia Tensor into cubic form.
  double p = - (Ixx + Iyy + Izz);
  double q =   Iyy*Izz + Iyy*Ixx + Izz*Ixx - (Ixy*Ixy + Iyz*Iyz + Ixz*Ixz);
  double r = - Ixx*Iyy*Izz + Ixx*Iyz*Iyz + Iyy*Ixz*Ixz +
    Izz*Ixy*Ixy - 2.*Ixy*Iyz*Ixz;

  // See CRC's Standard Mathematical Tables (19th edition), pp 103-105.
  // The substitution, y = x - p/3 converts  y^3 + p*y^2 + q*y + r = 0
  // to the form  x^3 + a*x + b = 0 .  Then, if b^2/4 + a^3/27 < 0 ,
  // there will be three real roots -- guaranteed since the Inertia Tensor
  // is symmetric.  A second substitution, x = m*cos(psi) , yields the roots.
  double a = (3.*q - p*p)/3.;
  double b = (2.*p*p*p - 9.*p*q + 27.*r)/27.;

  double rad_test = b*b/4. + a*a*a/27.;
  if ( (rad_test < 0.) && (Ixy != 0.) && (Ixz != 0.) && (Iyz != 0.) ) {
    // Update the weight and centroid
    m_weightSum  = weightSum;
    m_centroid   = centroid;
    m_centroid  /= weightSum;

    // Construct the roots, which are the principal moments.
    double m   = 2. * sqrt(-a/3.);
    double psi = acos( 3.*b/(a*m) ) / 3.;
    
    m_moment[0] = m * cos(psi) - p/3.;
    m_moment[1] = m * cos(psi + 2.*M_PI/3.) - p/3.;
    m_moment[2] = m * cos(psi + 4.*M_PI/3.) - p/3.;

    // Construct direction cosines; dircos for middle root is parallel to
    // longest principal axis.
    for(int iroot=0; iroot < 3; iroot++) 
      {
        double A = Iyz * (Ixx - m_moment[iroot]) - Ixy*Ixz;
        double B = Ixz * (Iyy - m_moment[iroot]) - Ixy*Iyz;
        double C = Ixy * (Izz - m_moment[iroot]) - Ixz*Iyz;
        double D = sqrt( 1. / ( 1./(A*A) + 1./(B*B) + 1./(C*C) ) ) / C;

        m_axis[iroot] = Vector(D*C/A, D*C/B, D);

        // Set axis to "point up"
        if ( m_axis[iroot].z() < 0. ) {
          m_axis[iroot] = -m_axis[iroot];
        }
      }

    // Calculate the covariance on the primary axis
    m_axisErr = calcCovariance(m_axis[0]);

    // Second loop to get the chisquare (residuals about principal axis, through centroid,
    // using input weight), the full cluster length, the skewness and the fraction of
    // energy inside the core cylinder.
    double xdir = m_axis[1].x();
    double ydir = m_axis[1].y();
    double absxdir = fabs(xdir);
    double absydir = fabs(ydir);
    chisq = 0.; 
    double tmin = 9999.;
    double tmax = -9999.;
    double tave = 0.;
    double tvar = 0.;
    double tskew = 0.;
    double coreEnergy = 0.;
    for(vecIter = dataVec.begin(); vecIter != dataVec.end(); vecIter++)
      {
        CalMomentsData& dataPoint = *vecIter;
        double distToAxis = dataPoint.calcDistToAxis(m_centroid, m_axis[1]);
        double weight = dataPoint.getWeight();
        double t = dataPoint.calcCoordAlongAxis(m_centroid, m_axis[1]);

        // This is a rough attempt to correct for the gaps in the CAL modules;
        // the gaps in the x and y directions are treated separately so the whole
        // thing is really incorrect and, if these quantities will ever turn out
        // to be useful we'll probably have to devise something smarter, here
        // (i. e. use a real propagator).
        // Also the gap between towers (45.6 mm) should not be hard-coded but
        // retrieved from the appropriate service.
        int tower = dataPoint.getTower();
        if ( absxdir > 0.05 ) {
          t -= 45.6*(tower % 4)/xdir; // (tower % 4) counts the gaps in the x direction.
        }
        if ( absydir > 0.05 ) {
          t -= 45.6*(tower / 4)/ydir; // (tower / 4) counts the gaps in the y direction.
        }
        // End of the "embarassing" part of the code (Luca Baldini, Dec. 26, 2010).
        
        chisq += dataPoint.getWeight()*distToAxis*distToAxis;
        if ( t < tmin ) {
          tmin = t;
        }
        if ( t > tmax ) {
          tmax = t;
        }
        tave  += weight * t;
        tvar  += weight * t*t;
        tskew += weight * t*t*t;
        if ( distToAxis < coreRadius ) {
          coreEnergy += dataPoint.getWeight();
        }
      }
    // Normalize all this garbage and refer to centroid.
    tave  /= m_weightSum;
    tvar  /= m_weightSum;
    tvar  -= tave*tave;
    tskew /= m_weightSum;
    tskew -= (3*tvar*tave + tave*tave*tave);
    if ( tvar > 0. ) {
      tskew /= pow(tvar, 1.5);
    }

    // Scale the chisquare by number of data points.
    chisq /= weightSum * dataVec.size();

    // Scale the distance between the two extreme points in the longitudinal
    // shower profile in order to convert the length from mm to X0.
    // The X0 value for the CsI(Tl) is taken from the PDG Review of Particle Physics (2008).
    m_fullLength = (tmax - tmin)/18.6;

    // Final calculations to return moment of principal axis and average of other two.
    // Note that the normalization to the sum of weights, which used to be done in the
    // AnalysisNtuple package, is done here (Luca Baldini, Dec 23, 2010).
    double longMag1  = fabs(m_moment[0]);
    double longMag2  = fabs(m_moment[2]);
    double transMag  = fabs(m_moment[1]);
    m_longRms        = sqrt((longMag1 + longMag2) / (2.*m_weightSum));
    m_transRms       = sqrt(transMag/m_weightSum);
    m_longRmsAsym    = (longMag1 - longMag2)/(longMag1 + longMag2);
    m_longSkewness   = tskew;
    m_coreEnergyFrac = coreEnergy/m_weightSum;
  }

  // This is when the radix test fails.
  else {
    chisq = -1.;
  }
  
  return chisq;
}

double CalMomentsAnalysis::doIterativeMomentsAnalysis(CalMomentsDataVec dataVec,
                                                      const Point& inputCentroid,
                                                      double transScaleFactor,
                                                      double transScaleFactorBoost,
                                                      double coreRadius)
{
  // First reset the class members keeping track of the iteration stats.
  m_numIterations    = 0;
  m_numDroppedPoints = 0;
  
  // Initialize some local variables.
  double chisq   = -1.;
  bool iterate   = true;
  Point centroid = inputCentroid;

  // Iterate until either failure (chisq < 0) or all remaining data points 
  // are within "range" of axis
  while(iterate)
    {
      // Do the standard moments analysis
      double localChisq = doMomentsAnalysis(dataVec, centroid, coreRadius);

      // Make sure it didn't fail on this iteration
      if ( localChisq < 0. ) {
        break;
      }

      // Update global chi-square to pick up this iteration's value
      chisq = localChisq;

      // Update the centroid for subsequent passes
      centroid = getCentroid();

      // Get the transverse moment
      double transRms = getTransRms();

      // Sort the data points by their distance from the principal axis
      std::sort(dataVec.begin(), dataVec.end());

      // Assume all data within range
      iterate = false;

      m_numIterations++;

      // Check the last element in the now sorted list of data points
      // and see if it is out of range
      while(!dataVec.empty())
        {
          CalMomentsData& momentsData = dataVec.back();
          
          // If out of range drop this point and loop back to check again
          if ( momentsData.getDistToAxis() > transScaleFactor * transRms ) {
            dataVec.pop_back();
            iterate = true;
            m_numDroppedPoints++;
          }
          else break;
        }

      // Make it harder to drop points on each iteration by boosting the scale factor
      transScaleFactor *= transScaleFactorBoost;

      // Finally, make sure have enough points remaining to proceed!
      if (dataVec.size() < 3) break;
    }

  return chisq;
}

CLHEP::HepMatrix CalMomentsAnalysis::calcCovariance(Vector momAxis) {
  // ADW: 11/21/11 Place holder for calculating the covariant errors from 
  // the moments analysis. For the time being, just use a
  // uniform (and naive) 2 degree error on the moment direction

  /// Hardcoded errors on the angles (to start with)
  double dtheta = 2*(M_PI/180.);
  double dphi = 2*(M_PI/180.);

  double theta = momAxis.theta();
  double phi = momAxis.phi();

  /// Define all the trig
  double st = sin(theta) ;  double st2 = st*st ; 
  double ct = cos(theta) ;  double ct2 = ct*ct ; 
  double sp = sin(phi)	 ;  double sp2 = sp*sp ; 
  double cp = cos(phi)	 ;  double cp2 = cp*cp ; 

  double dt2 = dtheta*dtheta ;
  double dp2 = dphi*dphi     ;

  CLHEP::HepMatrix cov(3,3,1);

  /// Set the approximate covariance
  cov(1, 1) = ct2*cp2*dt2 + st2*sp2*dp2	     ;
  cov(2, 2) = ct2*sp2*dt2 + st2*cp2*dp2	     ;
  cov(3, 3) = st2*dt2 		       	     ;
  cov(1, 2) = ct2*cp*sp*dt2 - st2*cp*sp*dp2  ;
  cov(1, 3) = -st*ct*cp*dt2	       	     ;
  cov(2, 3) = -st*ct*sp*dt2	       	     ;
  cov(2, 1) = cov(1, 2)		       	     ;
  cov(3, 1) = cov(1, 3)		       	     ;
  cov(3, 2) = cov(2, 3)		       	     ;

  /// Just a test to see everything is working
  //cov(1,1) = 1;
  //cov(2,2) = 2;
  //cov(3,3) = 3;

  return cov;
}
