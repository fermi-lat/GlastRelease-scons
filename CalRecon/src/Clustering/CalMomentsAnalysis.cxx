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
    // full calculation based on xtal position: 
    // update both centroid and direction matrix
    calcCovariance(dataVec, iniCentroid);
    // simple implementation based on theta/phi error
    //calcCovarianceAxisSimple(m_axis[1]); // REMEMBER: direction is Axis[1]
    //calcCovarianceCentroidSimple(dataVec, m_axis[1], iniCentroid); // Centroid part

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



void CalMomentsAnalysis::calcCovariance(CalMomentsDataVec& dataVec,
                                        const Point& centroid)
{
  // This is Bill McConville implementation of cov matrix calculation
  // starting from uncertainties in xtal position and error.
  // It follow the formalism in CalRecon/doc/moments_analysis.pdf 
  Vector mom=CalMomentsAnalysis::getMoments();
  double *L0=&mom[0];
  double *L1=&mom[1];
  double *L2=&mom[2];
  //
  //Define S, the transformation matrix that takes you from instrument coordinates to
  //the Primary axes coordinates.  This is constructed using the (transposed) eigenvectors obtained
  //from the moment analysis as the rows
  CLHEP::HepMatrix S(3,3,0);
  for(int i=0;i<3;i++)
  {
    for(int j=0;j<3;j++)
    {
      S(i+1,j+1) = CalMomentsAnalysis::getAxis(i)[j];
    }
  }
  //If S was defined correctly, you should retrieve the diagonalized matrix after tranforming on I
  //CLHEP::HepMatrix transform_I=S*CalMomentsAnalysis::getInertiaTensor()*S.T();
  
  //Define Sp, used to obtain the transformation matrix K
  CLHEP::HepMatrix Sp(9,3,0);
  Sp(1,2) = S(3,1);
  Sp(1,3) = -1*S(2,1);
  Sp(2,1) = -1*S(3,1);
  Sp(2,3) = S(1,1);
  Sp(3,1) = S(2,1);
  Sp(3,2) = -1*S(1,1);
  Sp(4,2) = S(3,2);
  Sp(4,3) = -1*S(2,2);
  Sp(5,1) = -1*S(3,2);
  Sp(5,3) = S(1,2);
  Sp(6,1) = S(2,2);
  Sp(6,2) = -1*S(1,2);
  Sp(7,2) = S(3,3);
  Sp(7,3) = -1*S(2,3);
  Sp(8,1) = -1*S(3,3);
  Sp(8,3) = S(1,3);
  Sp(9,1) = S(2,3);
  Sp(9,2) = -1*S(1,3);

  //Define D, used to rearrange the vec operator
  CLHEP::HepMatrix D(6,9,0);
  D(1,1) = 1;
  D(2,5) = 1;
  D(3,9) = 1;
  D(4,2) = 0.5;
  D(4,4) = 0.5;
  D(5,3) = 0.5;
  D(5,7) = 0.5;
  D(6,6) = 0.5;
  D(6,8) = 0.5;

  //Define Dplus, the pseudo inverse of D
  CLHEP::HepMatrix Dplus(9,6,0);
  Dplus(1,1) = 1;
  Dplus(2,4) = 1;
  Dplus(3,5) = 1;
  Dplus(4,4) = 1;
  Dplus(5,2) = 1;
  Dplus(6,6) = 1;
  Dplus(7,5) = 1;
  Dplus(8,6) = 1;
  Dplus(9,3) = 1;

  //Define Gplus, also used for finding K
  double g1,g2,g3;
  g1 = 0.5/(*L2 - *L1);
  g2 = 0.5/(*L0 - *L2);
  g3 = 0.5/(*L1 - *L0);
  CLHEP::HepMatrix Gplus(6,9,0);
  Gplus(1,1) = 1;
  Gplus(2,5) = 1;
  Gplus(3,9) = 1;
  Gplus(4,6) = g1;
  Gplus(4,8) = g1;
  Gplus(5,3) = g2;
  Gplus(5,7) = g2;
  Gplus(6,2) = g3;
  Gplus(6,4) = g3;

  //Get the Kronecker product of S and S, which has to be done by hand...
  CLHEP::HepMatrix Skron(9,9,0);
  for(int i=0; i<3; i++)
  {
    for(int j=0; j<3; j++)
    {
      for(int k=1; k<4; k++)
      {
        for(int l=1; l<4; l++)
	{
	  Skron(3*i+k,3*j+l)=S(i+1,j+1)*S(k,l);
	}
      }
    }
  }
  
  //Define Finverse, the final matrix needed to form K
  CLHEP::HepMatrix Finverse = Gplus * Skron * Dplus;

  CLHEP::HepMatrix K_left(12,6,0);
  for(int i=1;i<10;i++)
  {
    for(int j=1;j<4;j++)
    {
      K_left(i,j+3)=Sp(i,j);
    }
  }
  K_left(10,1)=1;
  K_left(11,2)=1;
  K_left(12,3)=1;

  //calculate K (finally)
  CLHEP::HepMatrix K=K_left*Finverse;

  
  //initialize covariance matrix components of inertia tensor
  double cIxx_xx = 0.0;
  double cIxx_yy = 0.0;
  double cIxx_zz = 0.0;
  double cIxx_xy = 0.0;
  double cIxx_xz = 0.0;
  double cIxx_yz = 0.0;
	
  double cIyy_yy = 0.0;
  double cIyy_zz = 0.0;
  double cIyy_xy = 0.0;
  double cIyy_xz = 0.0;
  double cIyy_yz = 0.0;

  double cIzz_zz = 0.0;
  double cIzz_xy = 0.0;
  double cIzz_xz = 0.0;
  double cIzz_yz = 0.0;	

  double cIxy_xy = 0.0;
  double cIxy_xz = 0.0;
  double cIxy_yz = 0.0;

  double cIxz_xz = 0.0;
  double cIxz_yz = 0.0;

  double cIyz_yz = 0.0;


  //initialize the coefficients for the energy centroid covariance matrix.
  //This is a diagonal matrix in instrument coordinates, so we only need
  //three terms
  double cCxx = 0.0;
  double cCyy = 0.0;
  double cCzz = 0.0;
  double wtot = 0.0; //set the total weight to zero

  // Loop through each of the hits and add the errors to the respective
  // components of the covariance matrix
  double x,y,z,w,dx,dy,dz,dw,x2,y2,z2,w2,dx2,dy2,dz2,dw2,ww;
  double tlen = 27; // from xtal geometry
  double llen = 326;
  double zlen = 20;
  double sqrt12 = sqrt(12.);
  CalMomentsDataVec::iterator vecIter;
  // cs: get wtot first
  for(vecIter = dataVec.begin(); vecIter != dataVec.end(); vecIter++) 
    { const CalMomentsData& dataPoint = *vecIter;
      wtot += dataPoint.getWeight(); 
    }
  // real loop for errors
  for(vecIter = dataVec.begin(); vecIter != dataVec.end(); vecIter++)
    {
      const CalMomentsData& dataPoint = *vecIter;

      Vector hit = dataPoint.getPosition() - centroid;
      x = hit.x();
      y = hit.y();
      z = hit.z();
      w = dataPoint.getWeight();
      x2 = x*x;
      y2 = y*y;
      z2 = z*z;
      w2 = w*w;  

      Vector main_axis = CalMomentsAnalysis::getAxis(1);      
      
      // Error values based on crystal dimensions (Atwood et al 2009). 
      /*
      double xlen, ylen, zlen;
      // Set the Xtal uncertainties according to WMC
      if (dataPoint.isx()){
	  xlen = 326;
	  ylen = 27;
	  zlen = 20;
	  if (ylen * fabs(main_axis[2]) > zlen * fabs(main_axis[1])){ 
	      if (xlen * fabs(main_axis[2]) > zlen * fabs(main_axis[0])){
	          dx = 5.0 + 0.5 * (zlen / fabs(main_axis[2])) * fabs(main_axis[0]);
	      }
	      else{
		  dx = xlen / 2;
	      }
	  }
	  else{
	      if (xlen * fabs(main_axis[1]) > ylen *fabs( main_axis[0])){
		  dx = 5.0 + 0.5 * (ylen / fabs(main_axis[1])) * fabs(main_axis[0]);
	      }
	      else{
		  dx = xlen / 2;
	      }
	  }
          dy = ylen / 2;
          dz = zlen / 2;
      }
      else{
	  xlen = 27;
	  ylen = 326;
	  zlen = 20;
	  if (xlen * fabs(main_axis[2]) > zlen * fabs(main_axis[0])){ 
	      if (ylen * fabs(main_axis[2]) > zlen * fabs(main_axis[1])){
	          dy = 5.0 + 0.5 * (zlen / fabs(main_axis[2])) * fabs(main_axis[1]);
	      }
	      else{
		  dy = ylen / 2;
	      }
	  }
	  else{
	      if (ylen * fabs(main_axis[0]) > xlen * fabs(main_axis[1])){
		  dy = 5.0 + 0.5 * (xlen / fabs(main_axis[0])) * fabs(main_axis[1]);
	      }
	      else{
		  dy = ylen / 2;
	      }
	  }
          dx = xlen / 2;
          dz = zlen / 2;
      }
      dw = 0.1*w; // ERROR ON XTAL ENERGY SET TO  10%
      */

      // CS: Xtal uncertainties based on a different parametrization
      if (dataPoint.isx()){
        // Longitudinal Position Error (from a quick chat with Eric G.)
        // 10 mm for a MIP (~11MeV) than decrease as 1/sqrt(E) up to 1 mm at high E than stable
        // param as sqrt([0]*[0]/w + [1]*[1]); with [0] = 31; and [1] = 0.8
        // Johan suggested 1 + 9/sqrt(w/10)
        dx =  1. + 9./sqrt(w/10.);
        dy = tlen / sqrt12;
      }
      else{
        dy = 1. + 9./sqrt(w/10.); 
        dx = tlen / sqrt12;
      }
      dz = zlen / sqrt12;
      dw = 0.01*w; // ERROR ON XTAL ENERGY SET TO 1%
   
      
      //
      dx2 = dx*dx;
      dy2 = dy*dy;
      dz2 = dz*dz;
      dw2 = dw*dw;
      ww  = (1. - w/wtot)*(1. - w/wtot);

      //Ixx-others
      
      cIxx_xx +=  4*w2 * (y2 * dy2 + z2 * dz2) + pow((y2 + z2),2)  * dw2;
      cIxx_yy +=  4*w2 * z2 * dz2         + (y2 + z2)*(x2 + z2) * dw2;
      cIxx_zz +=  4*w2 * y2 * dy2         + (y2 + z2)*(x2 + y2) * dw2;
      cIxx_xy += -2*w2 * x*y  * dy2         - (y2 + z2)*(x*y)        * dw2;
      cIxx_xz += -2*w2 * x*z  * dz2         - (y2 + z2)*(x*z)        * dw2;
      cIxx_yz += -2*w2 * y*z  * (dy2+dz2) - (y2 + z2)*(y*z)        * dw2;

      // Iyy-Others
      cIyy_yy +=  4*w2 * (x2 * dx2 + z2 * dz2) + pow((x2 + z2),2)   * dw2;
      cIyy_zz +=  4*w2 * x2 * dx2	      + (x2 + z2)*(x2 + y2) * dw2; 
      cIyy_xy += -2*w2 * x*y  * dx2	      - (x2 + z2)*(x*y)	    * dw2;	   
      cIyy_xz += -2*w2 * x*z  * (dx2+dz2) - (x2 + z2)*(x*z)	    * dw2;	    
      cIyy_yz += -2*w2 * y*z  * dz2	      - (x2 + z2)*(y*z)	    * dw2;

      // Izz-Others
      cIzz_zz += 4*w2 * (x2 * dx2 + y2 * dy2) + pow((x2 + y2),2)    * dw2;
      cIzz_xy += -2*w2 * x*y  * (dx2+dy2) - (x2 + y2)*(x*y)	    * dw2;
      cIzz_xz += -2*w2 * x*z  * dx2	      - (x2 + y2)*(x*z)	    * dw2;
      cIzz_yz += -2*w2 * y*z  * dy2	      - (x2 + y2)*(y*z)	    * dw2;

      // Ixy-Others
      cIxy_xy += w2 * (y2 * dx2 + x2 * dy2)   + x2 * y2         * dw2;
      cIxy_xz += w2 * y*z * dx2                     + x2 * y*z          * dw2;
      cIxy_yz += w2 * x*z * dy2                     + x * y2 * z        * dw2;

      // Ixz-Others
      cIxz_xz += w2 * (z2 * dx2 + x2 * dz2)   + x2 * z2         * dw2;
      cIxz_yz += w2 * x*y * dz2                     + x * y * z2        * dw2;
	    
      // Iyz-Iyz
      cIyz_yz += w2 * (y2 * dz2 + z2 * dy2)   + y2 * z2         * dw2;
      
      //add the errors in the energy centroid
      cCxx += w2 * dx2 + x2 * dw2 * ww;
      cCyy += w2 * dy2 + y2 * dw2 * ww;
      cCzz += w2 * dz2 + z2 * dw2 * ww;
      //cs: done in prev loop //wtot+=w;

    }
  //centroid error coefficients must be divided by total weight squared to get final answer
  cCxx/=pow(wtot,2);
  cCyy/=pow(wtot,2);
  cCzz/=pow(wtot,2);

  //create the 3x3 matrix for the error in the centroid
  CLHEP::HepMatrix centroidCovMatrix(3,3,0);
  centroidCovMatrix(1,1)=cCxx;
  centroidCovMatrix(2,2)=cCyy;
  centroidCovMatrix(3,3)=cCzz;

  CLHEP::HepMatrix VdICovMatrix(6,6,0);
  VdICovMatrix(1,1) = cIxx_xx;
  VdICovMatrix(1,2) = cIxx_yy;
  VdICovMatrix(1,3) = cIxx_zz;
  VdICovMatrix(1,4) = cIxx_xy;
  VdICovMatrix(1,5) = cIxx_xz;
  VdICovMatrix(1,6) = cIxx_yz;
  VdICovMatrix(2,1) = cIxx_yy;
  VdICovMatrix(2,2) = cIyy_yy;
  VdICovMatrix(2,3) = cIyy_zz;
  VdICovMatrix(2,4) = cIyy_xy;
  VdICovMatrix(2,5) = cIyy_xz;
  VdICovMatrix(2,6) = cIyy_yz;
  VdICovMatrix(3,1) = cIxx_zz;
  VdICovMatrix(3,2) = cIyy_zz;
  VdICovMatrix(3,3) = cIzz_zz;
  VdICovMatrix(3,4) = cIzz_xy;
  VdICovMatrix(3,5) = cIzz_xz;
  VdICovMatrix(3,6) = cIzz_yz;
  VdICovMatrix(4,1) = cIxx_xy;
  VdICovMatrix(4,2) = cIyy_xy;
  VdICovMatrix(4,3) = cIzz_xy;
  VdICovMatrix(4,4) = cIxy_xy;
  VdICovMatrix(4,5) = cIxy_xz;
  VdICovMatrix(4,6) = cIxy_yz;
  VdICovMatrix(5,1) = cIxx_xz;
  VdICovMatrix(5,2) = cIyy_xz;
  VdICovMatrix(5,3) = cIzz_xz;
  VdICovMatrix(5,4) = cIxy_xz;
  VdICovMatrix(5,5) = cIxz_xz;
  VdICovMatrix(5,6) = cIxz_yz;
  VdICovMatrix(6,1) = cIxx_yz;
  VdICovMatrix(6,2) = cIyy_yz;
  VdICovMatrix(6,3) = cIzz_yz;
  VdICovMatrix(6,4) = cIxy_yz;
  VdICovMatrix(6,5) = cIxz_yz;
  VdICovMatrix(6,6) = cIyz_yz;
  /*
  std::cout << "\nK Matrix:\n"<< K << "\n" << std::endl;
  std::cout << "\nVd Matrix:\n"<< VdICovMatrix << "\n" << std::endl;
  std::cout << "\nK.T Matrix:\n"<< K.T() << "\n" << std::endl;

  std::cout << "\nK*nVd Matrix:\n"<< K* VdICovMatrix << "\n" << std::endl;
  */

  CLHEP::HepMatrix errorCovarianceMatrix = K * VdICovMatrix * K.T();

  //this creates a covariance matrix for the principle moment
  CLHEP::HepMatrix cov(3,3,0);
  cov(1,1) = errorCovarianceMatrix(2,2);
  cov(1,2) = errorCovarianceMatrix(2,5);
  cov(1,3) = errorCovarianceMatrix(2,8);
  cov(2,1) = errorCovarianceMatrix(5,2);
  cov(2,2) = errorCovarianceMatrix(5,5);
  cov(2,3) = errorCovarianceMatrix(5,8);
  cov(3,1) = errorCovarianceMatrix(8,2);
  cov(3,2) = errorCovarianceMatrix(8,5);
  cov(3,3) = errorCovarianceMatrix(8,8);
  
  // Write global variable for error matrices, axis and centroid
  m_axisErr = cov;
  m_centroidErr = centroidCovMatrix;

}



void CalMomentsAnalysis::calcCovarianceAxisSimple(Vector momAxis) {
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
  double sp = sin(phi)         ;  double sp2 = sp*sp ; 
  double cp = cos(phi)         ;  double cp2 = cp*cp ; 

  double dt2 = dtheta*dtheta ;
  double dp2 = dphi*dphi     ;

  CLHEP::HepMatrix cov(3,3,1);

  /// Set the approximate covariance
  cov(1, 1) = ct2*cp2*dt2 + st2*sp2*dp2             ;
  cov(2, 2) = ct2*sp2*dt2 + st2*cp2*dp2             ;
  cov(3, 3) = st2*dt2                                     ;
  cov(1, 2) = ct2*cp*sp*dt2 - st2*cp*sp*dp2  ;
  cov(1, 3) = -st*ct*cp*dt2                            ;
  cov(2, 3) = -st*ct*sp*dt2                            ;
  cov(2, 1) = cov(1, 2)                                    ;
  cov(3, 1) = cov(1, 3)                                    ;
  cov(3, 2) = cov(2, 3)                                    ;

  /// Write the error matrix
  m_axisErr = cov;

  /// Just a test to see everything is working
  //cov(1,1) = 1;
  //cov(2,2) = 2;
  //cov(3,3) = 3;
  //return cov;

}

void CalMomentsAnalysis::calcCovarianceCentroidSimple(CalMomentsDataVec& dataVec, Vector momAxis,
                                                      const Point& centroid) {
  // CS This is a placeholder for Centroid error matrix calculation
  // from WMC implementation in which centroid cov matrix is diagonal
  
  //initialize the coefficients for the energy centroid covariance matrix.
  //This is a diagonal matrix in instrument coordinates, so we only need
  //three terms
  double cCxx = 0.0;
  double cCyy = 0.0;
  double cCzz = 0.0;
  double wtot = 0.0; //set the total weight to zero

  // Loop through each of the hits and add the errors to the respective
  // components of the covariance matrix
  double x,y,z,w,dx,dy,dz,dw,x2,y2,z2,w2,dx2,dy2,dz2,dw2;
  double tlen = 27;
  double llen = 326;
  double zlen = 20;
  double sqrt12 = sqrt(12.);
  CalMomentsDataVec::iterator vecIter;
  for(vecIter = dataVec.begin(); vecIter != dataVec.end(); vecIter++)
    {
      const CalMomentsData& dataPoint = *vecIter;
      
      Vector hit = dataPoint.getPosition() - centroid;
      x = hit.x();
      y = hit.y();
      z = hit.z();
      w = dataPoint.getWeight();
      x2 = x*x;
      y2 = y*y;
      z2 = z*z;
      w2 = w*w;  
      
      
      // Fill position errors
      // Error values are currently based off of crystal dimensions (Atwood et al 2009). 
 
      if (dataPoint.isx()){
        // Longitudinal Position Error (from a quick chat with Eric G.)
        // 10 mm for a MIP (~11MeV) than decrease as 1/sqrt(E) up to 1 mm at high E than stable
        // param as sqrt([0]*[0]/w + [1]*[1]); with [0] = 31; and [1] = 0.8
        // Johan suggested 1 + 9/sqrt(w/10)
        dx =  1. + 9./sqrt(w/10.);
        dy = tlen / sqrt12;
        dz = zlen / sqrt12;
      }
      else{
        dy = 1. + 9./sqrt(w/10.); //sqrt(31.*31./w + 0.8*0.8); 
        dx = tlen / sqrt12;
        dz = zlen / sqrt12;
      }
      dw = 0.01*w; // ERROR ON XTAL ENERGY SET TO 1%
      dx2 = dx*dx;
      dy2 = dy*dy;
      dz2 = dz*dz;
      dw2 = dw*dw;
      //add the errors in the energy centroid
      cCxx += w2 * dx2 + x2 * dw2;
      cCyy += w2 * dy2 + y2 * dw2;
      cCzz += w2 * dz2 + z2 * dw2;
      wtot+=w;

    }
  //centroid error coefficients must be divided by total weight squared to get final answer
  cCxx/=pow(wtot,2);
  cCyy/=pow(wtot,2);
  cCzz/=pow(wtot,2);
  
  //create the 3x3 matrix for the error in the centroid
  m_centroidErr(1,1)=cCxx;
  m_centroidErr(2,2)=cCyy;
  m_centroidErr(3,3)=cCzz;
}
