
#include "MomentsClusterInfo.h"
#include "src/Utilities/CalException.h"

#include "CalUtil/CalDefs.h"
#include "Event/Digi/CalDigi.h"

#include "idents/CalXtalId.h"

namespace 
{
    #ifdef WIN32
    #include <float.h> // used to check for NaN
    #else
    #include <cmath>
    #endif
     
    bool isFinite(double val) 
    {
        using namespace std; // should allow either std::isfinite or ::isfinite
        #ifdef WIN32
            return (_finite(val)!=0);  // Win32 call available in float.h
        #else
            return (isfinite(val)!=0); // gcc call available in math.h
        #endif
    }
} // anom namespace
 
Point MomentsClusterInfo::m_fit_p_layer[8];
Vector MomentsClusterInfo::m_fit_v_layer[8];
double MomentsClusterInfo::m_fit_e_layer[8];
double MomentsClusterInfo::m_fit_errp_layer[8];
double MomentsClusterInfo::m_fit_chisq = 0;
double MomentsClusterInfo::m_fit_xcentroid = 0;
double MomentsClusterInfo::m_fit_ycentroid = 0;
double MomentsClusterInfo::m_fit_zcentroid = 0;
double MomentsClusterInfo::m_fit_xdirection = 0;
double MomentsClusterInfo::m_fit_ydirection = 0;
double MomentsClusterInfo::m_fit_zdirection = 0;
int MomentsClusterInfo::m_fit_option = 0;

/// Constructor
MomentsClusterInfo::MomentsClusterInfo(ICalReconSvc* calReconSvc, ICalCalibSvc *calCalibSvc) :
  m_calReconSvc(calReconSvc),
  m_calCalibSvc(calCalibSvc),  
  m_p0(0.,0.,0.),
  m_p1(0.,0.,0.)
{
  m_calnLayers = m_calReconSvc->getCalNLayers();
  // Minuit object
  m_minuit = new TMinuit(5);
  //Sets the function to be minimized
  m_minuit->SetFCN(fcncal);
  return;
}

/// This makes CalClusters out of associated CalXtalRecData pointers
Event::CalCluster* MomentsClusterInfo::fillClusterInfo(const XtalDataList* xTalVec)
{
  // Create an output cluster.
  Event::CalCluster* cluster = 0;
  
  if ( !xTalVec->empty() ) {
    cluster = new Event::CalCluster(0);
    cluster->setProducerName("MomentsClusterInfo");
    cluster->clear();
    
    // Do Philippe's fit.
    fitDirectionCentroid(xTalVec);

    // Perform longitudinal position and energy corrections for saturated crystals
    CorrectSaturatedXtalLongPositionAndEnergy(xTalVec);

    // Calculate basic cluster properties.
    fillLayerData(xTalVec, cluster);
    
    // Do the iterative moments analysis.
    fillMomentsData(xTalVec, cluster);
  }
  return cluster;
}


void MomentsClusterInfo::fillLayerData(const XtalDataList* xTalVec,
                                       Event::CalCluster* cluster)
{
  int numXtals           = xTalVec->size();
  double xtalRawEneSum   = 0.;
  double xtalCorrEneSum  = 0.;
  double xtalBadEneSum   = 0.;
  Point centroid(0., 0., 0.);
  std::vector<double> layerEnergy(m_calnLayers, 0.);
  std::vector<Vector> layerAve(m_calnLayers);
  std::vector<Vector> layerRms(m_calnLayers);

  // Use the fit centroid as the start point.
  m_p0 = Point(0.,0.,0.);

  // Loop over all crystals in the current cluster and compute barycenter
  // and various moments.
  XtalDataList::const_iterator xTalIter;
  
  for ( xTalIter = xTalVec->begin(); xTalIter != xTalVec->end(); xTalIter++ )
    {
      // Get pointer to the reconstructed data for given crystal.
      Event::CalXtalRecData* recData = *xTalIter;
      // Get reconstructed values.
      double xtalEne = recData->getEnergy();
      Vector xtalPos = recData->getPosition();
      int    layer   = (recData->getPackedId()).getLayer();
      // Check if bad position measurement.
      if ( !recData->isPositionGood() ) {
        xtalBadEneSum += xtalEne;
        continue;
      }

      // Offset the position to be closer to actual cluster center.
      xtalPos -= m_p0;
      
      // Update energy of corresponding layer
      layerEnergy[layer] += xtalEne;
      
      // Update average position of corresponding layer.
      Vector wpos = xtalEne*xtalPos;
      layerAve[layer] += wpos;
  
      // Vector containing squared coordinates, weighted by crystal energy.
      Vector wpos2(wpos.x()*xtalPos.x(), wpos.y()*xtalPos.y(),
                   wpos.z()*xtalPos.z());
      
      // Update quadratic spread, which is proportional to xtalEne;
      // this means, that position error in one crystal
      // is assumed to be 1/sqrt(xtalEne)
      layerRms[layer] += wpos2;
      
      // Update energy sum
      xtalRawEneSum += xtalEne;
      
      // update cluster position
      centroid += wpos;
    }
  
  // If energy sum is not zero normalize cluster position...
  if ( xtalRawEneSum > 0. ) {
    centroid /= xtalRawEneSum;
    centroid += m_p0;
    cluster->setStatusBit(Event::CalCluster::CENTROID);
  }
  // ...otherwise set cluster position to non-physical value
  else {
    centroid = Point(-1000., -1000., -1000.);
  }
    
  // Loop over calorimeter layers.
  for ( int i = 0; i < m_calnLayers; i++ )
    {
      // If energy in the layer is not zero finalize calculations...
      if( layerEnergy[i] > 0 ) {
        // Normalize position in the layer.
        layerAve[i] /= layerEnergy[i]; 
        // Normalize quadratic spread in the layer.
        layerRms[i] /= layerEnergy[i];

        // Vector containing the squared average position in each component
        Vector sqrLayer(layerAve[i].x()*layerAve[i].x(),
                        layerAve[i].y()*layerAve[i].y(),
                        layerAve[i].z()*layerAve[i].z());
        
        // the precision of transverse coordinate measurement
        // if there is no fluctuations: 1/sqrt(12) of crystal width
        Vector d;
        double csIWidth = m_calReconSvc->getCalCsIWidth();
        if ( i%2 == 1 ){
          d = Vector(csIWidth*csIWidth/12.,0.,0.);
        }
        else {
          d = Vector(0.,csIWidth*csIWidth/12.,0.);
        }
  
        // Subtracting the  squared average position and adding the square of
        // crystal width, divided by 12.
        layerRms[i] += d - sqrLayer;
        
        // Reset layerAve to detector coordinates
        layerAve[i] += m_p0;
      }
        
      // Otherwise reset position and spread Vectors.
      else {
        layerAve[i] = m_p0;
        layerRms[i] = m_p0;
      }

      // Fill the cluster layer data.
      Point layerPos(layerAve[i].x(), layerAve[i].y(), layerAve[i].z());
      
      // Set the data for the vector.
      Event::CalClusterLayerData layerData(layerEnergy[i], layerPos,
                                           layerRms[i]);

      cluster->push_back(layerData);
    }

  // Total energy in cluster includes energy in bad position crystals
  xtalRawEneSum += xtalBadEneSum;

  // Loop over all crystals in the current cluster again to compute the number
  // of truncated xtals, the number of saturated xtals and the moments of the
  // xtal energy distribution.
  int numTruncXtals      = 0;
  int numEneMomXtals     = 0;
  int numSaturatedXtals  = 0;
  double xtalEneMax      = 0.;
  double xtalEneMomSum   = 0.;
  double xtalEneMomSum2  = 0.;
  double xtalEneMomSum3  = 0.;
  double xtalsTruncFrac  = m_calReconSvc->getMciXtalsTruncFrac();
  double eneMomTruncFrac = m_calReconSvc->getMciEneMomTruncFrac();
  for( xTalIter = xTalVec->begin(); xTalIter != xTalVec->end(); xTalIter++ )
    {
      Event::CalXtalRecData* recData = *xTalIter;
      double xtalEne = recData->getEnergy();
      if ( xtalEne > xtalEneMax ) {
        xtalEneMax = xtalEne;
      }
      if ( xtalEne > (xtalsTruncFrac*xtalRawEneSum) ) {
        numTruncXtals ++;
      }
      if ( xtalEne > (eneMomTruncFrac*xtalRawEneSum) ) {
        numEneMomXtals ++;
        xtalEneMomSum  += xtalEne;
        xtalEneMomSum2 += xtalEne*xtalEne;
        xtalEneMomSum3 += xtalEne*xtalEne*xtalEne;
      }
      if ( recData->saturated() ) {
        numSaturatedXtals += 1;
      }
    }

  // Normalize the moments of the xtal energy distribution, if needed.
  double xtalEneRms = -1.;
  double xtalEneSkewness = 0.;
    
  if ( numEneMomXtals > 1 ) {
    xtalEneMomSum  /= numEneMomXtals;
    xtalEneMomSum2 /= numEneMomXtals;
    xtalEneMomSum3 /= numEneMomXtals;
    
    if ( xtalEneMomSum2 > 0. ) {
      xtalEneRms = xtalEneMomSum2 - xtalEneMomSum*xtalEneMomSum;
      xtalEneSkewness = (xtalEneMomSum3 - 3*xtalEneMomSum*xtalEneRms - 
                         xtalEneMomSum*xtalEneMomSum*xtalEneMomSum);
      if(xtalEneRms<=0)
        {
          xtalEneRms = 0;
          xtalEneSkewness = 0;
        }
      else
        {
          xtalEneRms = sqrt(xtalEneRms);
          xtalEneSkewness /= (xtalEneRms*xtalEneRms*xtalEneRms);
        }
    }
    if ( !isFinite(xtalEneRms) ) {
      throw CalException("MomentsClusterInfo: infinite xtalEneRms");
    }
    if ( !isFinite(xtalEneSkewness) ) {
      throw CalException("MomentsClusterInfo: infinite xtalEneSkewness");
    }
  }

  // Correct for Zero Supression using Truncated Xtal count.
  double zeroSupprEnergy = m_calReconSvc->getMciZeroSupprEnergy();
  xtalCorrEneSum = xtalRawEneSum + (zeroSupprEnergy*numTruncXtals);
  
  // Set the cluster data members: first the CalXtalsParams container...
  Event::CalXtalsParams xtalsParams(numXtals, numTruncXtals, numSaturatedXtals,
                                    xtalRawEneSum, xtalCorrEneSum, xtalEneMax,
                                    xtalEneRms, xtalEneSkewness, centroid);
  cluster->setXtalsParams(xtalsParams);

  // ...then we do create a minimal CalMomParams object in order to store the
  // centroid.
  // This is actually *not* the centroid of the moments analysis, but rather
  // the centroid calculated in the loop over the layers. The CalMomParams
  // member of the cluster object will be overwritten in the fillMomentsData()
  // methods. I don't think this is correct. When the moments analysis does
  // run all the way to the end it does not make any difference, as this
  // member is overwritten, but whan that does not happen we end up with a bogus
  // class member faking something different from what it really is.
  // I tested with a few thousands gammas and there's indeed a difference if I
  // get rid of this, meaning that the CalMomParams information is used without
  // checking the CalCluster status bit.
  // Even stranger, it also makes a difference for the fit direction/centroid,
  // thought for a tiny subset of the events.
  // I'll leave things as they are for the time being, but we might want to
  // revise this.
  //
  // Luca Baldini, January 18, 2011.
  Event::CalMomParams momParams(xtalCorrEneSum, 10*xtalCorrEneSum, centroid);
  cluster->setMomParams(momParams);
  
  // ...and eventually we set the fit parameters for the cluster. The reason
  // why we're doing it here is mainly historical, I guess, and it has to do
  // with the fact that the CalFitParams container has members for the energy
  // and the energy errors, though they're not calculated in the fit process.
  // We might want to refactor the code, at some point. We create the container
  // here because we might want to use the fit information in the moments
  // analysis.
  //
  // Also: it is important to check that the fit routine actually did
  // something by requiring that the number of fit layers is greater than 0,
  // since the fit variables are not properly reset and we carry memory of the
  // previous event if we don't do so. I'm not going to change Philippe's code,
  // just keep the check (modified) here.
  //
  // Luca Baldini, January 11, 2011.
  Event::CalFitParams fitParams(xtalCorrEneSum, 10*xtalCorrEneSum, centroid,
                                m_fit_nlayers, 0.);
  if ( m_fit_nlayers > 0 ) {
    fitParams = Event::CalFitParams(xtalCorrEneSum, 10*xtalCorrEneSum,
                                    m_fit_xcentroid, m_fit_ycentroid,
                                    m_fit_zcentroid,
                                    m_fit_xdirection, m_fit_ydirection,
                                    m_fit_zdirection,
                                    m_fit_nlayers, m_fit_chisq);
  }
  cluster->setFitParams(fitParams);
}

void MomentsClusterInfo::fillMomentsData(const XtalDataList* xtalVec,
                                         Event::CalCluster* cluster)
{
  // Grab the necessary xtalsParams...
  double xtalsCorrEnergySum = cluster->getXtalsParams().getXtalCorrEneSum();
  Point  xtalsCentroid = cluster->getXtalsParams().getCentroid();

  // ...and the fit params.
  Point fitCentroid = cluster->getFitParams().getCentroid();

  // Initialize empty centroid/axis to store the output of the moments analysis.
  Point momCentroid = Point(0.,0.,0.);
  Vector momAxis    = Vector(0.,0.,0.);

  // Create an empty vector of CalMomentsData object.
  CalMomentsDataVec momDataVec;
  momDataVec.clear();

  // Loop through the xtals setting the hits to analyze.
  XtalDataList::const_iterator xtalMax = xtalVec->end();
  XtalDataList::const_iterator xtalIter;

  for(xtalIter = xtalVec->begin(); xtalIter != xtalVec->end(); xtalIter++)
    {
      if ( xtalIter == xtalMax ) continue;
      
      CalMomentsData momData(*xtalIter);

      // This is done in MomentsClusterInfo::CorrectSaturatedXtalLongPositionAndEnergy (Ph. Bruel, 2012/12/02)
//       momData.applyFitCorrection(cluster->getFitParams(), m_calReconSvc);
//       // If the fit is reasonable, we might take advantage of it.
//       if ( (m_fit_nlayers >= 4) && (m_fit_zdirection > 0.) ) {
//         // Check whether the xtal is saturated.
//         if ( momData.saturated() ) {
//           momData.forceFitCorrection();
//         }
//         // If the longitudinal position is right on the edge of the xtal,
//         // or if the fit position is close to the xtal edge, use the fit
//         // position.
//         //if ( momData.checkStatusBit(CalMomentsData::LONG_POS_INVALID) ||
//         //     momData.checkStatusBit(CalMomentsData::FIT_POS_NEAR_EDGE) ) {
//         //  momData.enableFitCorrection();
//         //}
//       }    
      
      // Put the object into the vector.
      momDataVec.push_back(momData);
    }

  // Do the actual moments analysis.
  CalMomentsAnalysis momentsAnalysis;

  int numXtals = momDataVec.size();
  double transScaleFac = m_calReconSvc->getMaTransScaleFactor();
  double transScaleFacBoost = m_calReconSvc->getMaTransScaleFactorBoost();
  double coreRadius = m_calReconSvc->getMaCoreRadius();
  double chiSq = momentsAnalysis.doIterativeMomentsAnalysis(momDataVec,
                                                            xtalsCentroid,
                                                            transScaleFac,
                                                            transScaleFacBoost,
                                                            coreRadius);
  if ( chiSq >= 0. ) {
    // Get all the variables that need to be grabbed before the last iteration
    // (with all the xtals) is performed. This include the statistics on the
    // iterations, the centroid/axis, the longitudinal skewness and the full
    // cluster length.
    momCentroid = momentsAnalysis.getCentroid();
    momAxis = momentsAnalysis.getAxis();
    int numIterations = momentsAnalysis.getNumIterations();
    int numCoreXtals = numXtals - momentsAnalysis.getNumDroppedPoints();
    double longSkewness = momentsAnalysis.getLongSkewness();
    double fullLength = momentsAnalysis.getFullLength();

    // That done, recalculate the moments going back to using all the data
    // points but with the iterated moments centroid.
    if ( numIterations > 1 ) {
      chiSq = momentsAnalysis.doMomentsAnalysis(momDataVec, momCentroid,
                                                coreRadius);
    }
    
    // Extract the values for the moments with all hits present.
    double longRms  = momentsAnalysis.getLongRms();
    double transRms = momentsAnalysis.getTransRms();
    double longRmsAsym = momentsAnalysis.getLongRmsAsym();
    double coreEnergyFrac = momentsAnalysis.getCoreEnergyFrac();
    
    if ( !isFinite(longRms) ) {
      throw CalException("CalMomentsAnalysis: infinite longRms");
    }
    if ( !isFinite(transRms) ) {
      throw CalException("CalMomentsAnalysis: infinite transRms");
    }
    if ( !isFinite(longRmsAsym) ) {
      throw CalException("CalMomentsAnalysis: infinite longRmsAsym");
    }
    if ( !isFinite(longSkewness) ) {
      throw CalException("CalMomentsAnalysis: infinite longSkewness");
    }

    // Place holder for centroid error and axis covariance
    // ADW: Sept. 7, 2011
    CLHEP::HepMatrix momCentroidErr = momentsAnalysis.getCentroidErr();
    CLHEP::HepMatrix momAxisErr = momentsAnalysis.getAxisErr();

    // Store the information in the actual cluster: first the CalMomParams
    // container...
    //CLHEP::HepMatrix I_3_3(3, 3, 1);
    Event::CalMomParams momParams (xtalsCorrEnergySum, 10*xtalsCorrEnergySum,
                                   momCentroid, momCentroidErr, momAxis, momAxisErr,
                                   numIterations, numCoreXtals, numXtals,
                                   transRms, longRms, longRmsAsym, longSkewness,
                                   coreEnergyFrac, fullLength, -1., -1.);
    //momParams.fillStream(std::cout);
    cluster->setMomParams(momParams);

    // Set the relevant status bit and we're done!
    cluster->setStatusBit(Event::CalCluster::MOMENTS);
  }
  return;
}

int MomentsClusterInfo::getCentroidTransverseInfoOnly(const XtalDataList* xTalVec)
{
    //
    // Purpose and Method:
    // Find the centroid using only the transverse position information
    // This is needed when there are saturated crystals : we need the centroid and direction
    // in order to set the longitudinal position for saturated crystals before the moment analysis is run

    double mycentroid[3];
    double myenergy[3];
    int i;
    for(i=0;i<3;++i)
    {
        mycentroid[i] = 0;
        myenergy[i] = 0;
    }

    XtalDataList::const_iterator xTalIter;
    for(xTalIter = xTalVec->begin(); xTalIter != xTalVec->end(); xTalIter++)
    {
        // get pointer to the reconstructed data for given crystal
        Event::CalXtalRecData* recData = *xTalIter;
        double xtalEne = recData->getEnergy();                  // crystal energy
        int layer   = (recData->getPackedId()).getLayer();   // layer number
        Vector xtalPos   = recData->getPosition();         // Vector of crystal position
        if(layer%2==0)
        {
            myenergy[1] += xtalEne;
            mycentroid[1] += xtalEne*xtalPos.y();
        }
        else
        {
            myenergy[0] += xtalEne;
            mycentroid[0] += xtalEne*xtalPos.x();
        }
        myenergy[2] += xtalEne;
        mycentroid[2] += xtalEne*xtalPos.z();
    }
  
    if(myenergy[0]<=0 || myenergy[1]<=0) return 1; // we must have x AND y information

    for(i=0;i<3;++i)
        mycentroid[i] /= myenergy[i];

    m_p1 = Point(mycentroid[0],mycentroid[1],mycentroid[2]);

    return 0;
}

double MomentsClusterInfo::getsqdistancebetweenlines(Point p0, Vector v0, Point p1, Vector v1)
{
    Vector p0p1 = p0-p1;
    Point p2;
    double lambda = 1-(v0*v1)*(v0*v1);
    double mydist2 = 0;
    if(lambda==0) // xtal axis and main axis are parallel
    {
        mydist2 = p0p1*p0p1 - (p0p1*v1)*(p0p1*v1);
    }
    else
    {
        lambda = (-(p0p1*v0)+(p0p1*v1)*(v0*v1))/lambda;
        p2 = p0 + lambda*v0;
        p0p1 = p2-p1;
        mydist2 = p0p1*p0p1 - (p0p1*v1)*(p0p1*v1);
    }
    return mydist2;
}

double MomentsClusterInfo::compute_chi2_cal(double *par)
{
    if(m_fit_option>0)
    {
        m_fit_xcentroid = par[2];
        m_fit_ycentroid = par[3];
    }

    Point m_fit_p = Point(m_fit_xcentroid,m_fit_ycentroid,m_fit_zcentroid);
    m_fit_xdirection = sqrt(1-par[0]*par[0])*cos(par[1]);
    m_fit_ydirection = sqrt(1-par[0]*par[0])*sin(par[1]);
    m_fit_zdirection = par[0];
    Vector m_fit_v = Vector(m_fit_xdirection,m_fit_ydirection,m_fit_zdirection);

    m_fit_chisq = 0;
    int i;
    double mydist2 = 0;
    double xpos,ypos,zpos;
    for(i=0;i<8;++i)
    {
        if(m_fit_e_layer[i]<=0) continue;
        zpos = m_fit_p_layer[i].z();
        xpos = m_fit_p.x()+(zpos-m_fit_p.z())/m_fit_v.z()*m_fit_v.x();
        ypos = m_fit_p.y()+(zpos-m_fit_p.z())/m_fit_v.z()*m_fit_v.y();
        mydist2 = getsqdistancebetweenlines(m_fit_p_layer[i],m_fit_v_layer[i],m_fit_p,m_fit_v);
        m_fit_chisq += mydist2/64*m_fit_errp_layer[i];
    }

    return m_fit_chisq;
}

void MomentsClusterInfo::fcncal(int & , double *, double &f, double *par, int )
{
    f = compute_chi2_cal(par);
}


int MomentsClusterInfo::fitDirectionCentroid(const XtalDataList* xTalVec)
{
    //
    // Purpose and Method:

//     static Point m_fit_p_layer[8];
//     static Vector m_fit_v_layer[8];
//     static double m_fit_e_layer[8];
//     static double m_fit_chisq = 0;
//     static Point m_fit_p0;
//     static Point m_fit_p1;
//     static Point m_fit_v;

    m_fit_nlayers = 0;
    double totenergy = 0;
    double mypos[8][3];
    double mypos2[8][3];
    int i,j;
    for(i=0;i<8;++i)
    {
        m_fit_e_layer[i] = 0;
        for(j=0;j<3;++j)
        {
            mypos[i][j] = 0;
            mypos2[i][j] = 0;
        }
    }

    double mycentroid[3];
    double myenergy[3];
    for(i=0;i<3;++i)
    {
        mycentroid[i] = 0;
        myenergy[i] = 0;
    }

    XtalDataList::const_iterator xTalIter;
    for(xTalIter = xTalVec->begin(); xTalIter != xTalVec->end(); xTalIter++)
    {
        // get pointer to the reconstructed data for given crystal
        Event::CalXtalRecData* recData = *xTalIter;
        double xtalEne = recData->getEnergy();                  // crystal energy
        int layer   = (recData->getPackedId()).getLayer();   // layer number
        Vector xtalPos   = recData->getPosition();         // Vector of crystal position
        //
        totenergy += xtalEne;
        m_fit_e_layer[layer] += xtalEne;
        mypos[layer][0] += xtalEne*xtalPos.x();
        mypos[layer][1] += xtalEne*xtalPos.y();
        mypos[layer][2] += xtalEne*xtalPos.z();
        mypos2[layer][0] += xtalEne*xtalPos.x()*xtalPos.x();
        mypos2[layer][1] += xtalEne*xtalPos.y()*xtalPos.y();
        mypos2[layer][2] += xtalEne*xtalPos.z()*xtalPos.z();
        //
        if(layer%2==0)
        {
            myenergy[1] += xtalEne;
            mycentroid[1] += xtalEne*xtalPos.y();
        }
        else
        {
            myenergy[0] += xtalEne;
            mycentroid[0] += xtalEne*xtalPos.x();
        }
        myenergy[2] += xtalEne;
        mycentroid[2] += xtalEne*xtalPos.z();
    }

    if(myenergy[0]<=0 || myenergy[1]<=0) return 1; // we must have x AND y information
  
    for(i=0;i<3;++i)
        mycentroid[i] /= myenergy[i];

    m_fit_xcentroid = mycentroid[0];
    m_fit_ycentroid = mycentroid[1];
    m_fit_zcentroid = mycentroid[2];

    //
    m_fit_nlayers = 0;
    for(i=0;i<8;++i)
    {
        if(m_fit_e_layer[i]<=0) continue;
        ++m_fit_nlayers;
        for(j=0;j<3;++j)
        {
            mypos[i][j] /= m_fit_e_layer[i];
            mypos2[i][j] /= m_fit_e_layer[i];
            mypos2[i][j] -= mypos[i][j]*mypos[i][j];
            if(mypos2[i][j]<=0) mypos2[i][j] = 0;
        }
    }

    // need at least 4 layers
    if(m_fit_nlayers<4) return 1;

    for(i=0;i<8;++i)
    {
        if(m_fit_e_layer[i]<=0) continue;
        if(i%2==0)
        {
            m_fit_p_layer[i] = Point(0,mypos[i][1],mypos[i][2]);
            m_fit_v_layer[i] = Vector(1,0,0);
        }
        else
        {
            m_fit_p_layer[i] = Point(mypos[i][0],0,mypos[i][2]);
            m_fit_v_layer[i] = Vector(0,1,0);
        }
        m_fit_errp_layer[i] = m_fit_e_layer[i]/totenergy*(double)m_fit_nlayers;
    }

    //
    m_fit_option = 0;

    double arglist[10];
    int ierflg = 0;
    double amin,edm,errdef;
    int nvpar,nparx,icstat;
  
    // Set no output because of tuple output
    arglist[0] = -1;
    m_minuit->mnexcm("SET PRI", arglist ,1,ierflg);
  
    // idem with warnings ( minuit prints warnings
    // when the Hessian matrix is not positive )
    m_minuit->mnexcm("SET NOW", arglist ,1,ierflg);
  
    arglist[0] = 1;
    m_minuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
    // defines the strategy used by minuit
    // 1 is standard
    arglist[0] = 2;
    m_minuit->mnexcm("SET STR", arglist ,1,ierflg);        
  
    double vstart[4];
    double vstart2[4];
    double vstep[4];
  
    double par0,par1,par2,par3;
    double epar0,epar1,epar2,epar3;

    int nstep0 = 10;
    double par0min = 0;
    double par0max = 1;
    double par0step = (par0max-par0min)/(double)nstep0;
    int nstep1 = 4*nstep0;
    double par1min = -TMath::Pi();
    double par1max = TMath::Pi();
    double par1step = (par1max-par1min)/(double)nstep1;

    int imin = -1;
    int jmin = -1;
    double mymin = 99999999;
    double myval = 0;
    for(i=0;i<nstep0;++i)
    {
        for(j=0;j<nstep1;++j)
        {
            vstart2[0] = par0min+(0.5+(double)i)*par0step;
            vstart2[1] = par1min+(0.5+(double)j)*par1step;
            myval = compute_chi2_cal(vstart2);

            if(myval<mymin)
            {
                mymin = myval;
                imin = i;
                jmin = j;
            }
        }
    }

    vstart2[0] = par0min+(0.5+(double)imin)*par0step;
    vstart2[1] = par1min+(0.5+(double)jmin)*par1step;
    vstart[0] = vstart2[0];
    vstart[1] = vstart2[1];
    
    vstep[0] = par0step/2;
    vstep[1] = par1step/2;
  
    m_minuit->mnparm(0, "a1", vstart[0], vstep[0], par0min, par0max, ierflg);
    m_minuit->mnparm(1, "a2", vstart[1], vstep[1], par1min, par1max, ierflg);
  
    // Calls Migrad with 500 iterations maximum
    arglist[0] = 500;
    arglist[1] = 1.;
    m_minuit->mnexcm("MIGRAD", arglist ,2,ierflg); // fit with only two parameters : positions in x and y in the place z = mean z (direction fixed)

    m_minuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

    m_minuit->GetParameter(0,par0,epar0);
    m_minuit->GetParameter(1,par1,epar1);

    //
    m_fit_option = 1;

    // Clear minuit
    arglist[0] = 1;
    m_minuit->mnexcm("CLEAR", arglist ,1,ierflg);

    // Set no output because of tuple output
    arglist[0] = -1;
    m_minuit->mnexcm("SET PRI", arglist ,1,ierflg);
  
    // idem with warnings ( minuit prints warnings
    // when the Hessian matrix is not positive )
    m_minuit->mnexcm("SET NOW", arglist ,1,ierflg);
  
    arglist[0] = 1;
    m_minuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
    // defines the strategy used by minuit
    // 1 is standard
    arglist[0] = 2;
    m_minuit->mnexcm("SET STR", arglist ,1,ierflg);

    vstart[0] = par0;
    vstart[1] = par1;  
    vstart[2] = mycentroid[0];
    vstart[3] = mycentroid[1];
  
    vstep[0] = par0step/10;
    vstep[1] = par1step/10;
    vstep[2] = 5;
    vstep[3] = 5;
  
    m_minuit->mnparm(0, "a1", vstart[0], vstep[0], par0min, par0max, ierflg);
    m_minuit->mnparm(1, "a2", vstart[1], vstep[1], par1min, par1max, ierflg);
    m_minuit->mnparm(2, "a3", vstart[2], vstep[2], mycentroid[0]-100, mycentroid[0]+100, ierflg);
    m_minuit->mnparm(3, "a4", vstart[3], vstep[3], mycentroid[1]-100, mycentroid[1]+100, ierflg);
  
    // Calls Migrad with 500 iterations maximum
    arglist[0] = 500;
    arglist[1] = 1.;
    m_minuit->mnexcm("MIGRAD", arglist ,2,ierflg);  // fit with 4 parameters : theta and phi and positions in x and y

    m_minuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

    m_minuit->GetParameter(0,par0,epar0);
    m_minuit->GetParameter(1,par1,epar1);
    m_minuit->GetParameter(2,par2,epar2);
    m_minuit->GetParameter(3,par3,epar3);

    vstart2[0] = par0;
    vstart2[1] = par1;
    vstart2[2] = par2;
    vstart2[3] = par3;
    // Don't remove : the following line fills m_fit_xycentroid, m_fit_xyzdirection and m_fit_chisq !!!!!!!!!!!!!!!!!
    compute_chi2_cal(vstart2);

    // Clear minuit
    arglist[0] = 1;
    m_minuit->mnexcm("CLEAR", arglist ,1,ierflg);

    return 0;
}

int MomentsClusterInfo::CorrectSaturatedXtalLongPositionAndEnergy(const XtalDataList* xTalVec)
{
  if(!m_calCalibSvc) return 0;
  if( m_fit_nlayers==0 || m_fit_zdirection==0.) return 0;

  Point fitCentroid(m_fit_xcentroid,m_fit_ycentroid,m_fit_zcentroid);
  Vector fitAxis(m_fit_xdirection,m_fit_ydirection,m_fit_zdirection);

  double xtalcorposition;
  XtalDataList::const_iterator xTalIter;
  for(xTalIter = xTalVec->begin(); xTalIter != xTalVec->end(); xTalIter++)
    {
      // get pointer to the reconstructed data for given crystal
      Event::CalXtalRecData* recData = *xTalIter;
      if(!recData->saturated()) continue;
      if(recData->longPosCorrected() && recData->energyCorrected()) continue;
      //
      xtalcorposition = SetXtalLongPositionFromFitAxis(recData,fitCentroid,fitAxis);
      SetXtalEnergyFromPosition(recData,xtalcorposition);
    }
  return 0;
}

double MomentsClusterInfo::SetXtalLongPositionFromFitAxis(Event::CalXtalRecData* recData, Point fitCentroid, Vector fitAxis)
{
  int tower  = recData->getTower();
  int layer  = recData->getLayer();
  int column = recData->getColumn();
  // Otherwise carry on bravely and store the relevant information.
  double towerPitch = m_calReconSvc->getCaltowerPitch();
  double csILength  = m_calReconSvc->getCalCsILength();
  bool flightGeom   = m_calReconSvc->getCalFlightGeom();
  //
  double x = (recData->getPosition()).x();
  double y = (recData->getPosition()).y();
  double z = (recData->getPosition()).z();
  int towerIdy = tower / 4;
  int towerIdx = tower - 4*towerIdy;
  //
  double minCoord;
  double maxCoord;
  double distToEdge;
  // Xtals along the x coordinate first...
  double nocorposition = 0;
  double corposition = 0;
  double xtalcorposition = 0;
  if ( layer%2==0 ) 
    {
      nocorposition = x;
      minCoord = -1.5*towerPitch + towerPitch*(double)towerIdx - csILength/2;
      maxCoord = -1.5*towerPitch + towerPitch*(double)towerIdx + csILength/2;
      distToEdge = std::min( fabs(x - minCoord), fabs(x - maxCoord) );

      x = fitCentroid.x() + (z - fitCentroid.z()) / fitAxis.z() * fitAxis.x();
      distToEdge = std::min( fabs(x - minCoord), fabs(x - maxCoord) );

      if ( x < minCoord ) 
        {
          x = minCoord;
        }
      if ( x > maxCoord ) 
        {
          x = maxCoord;
        }
      corposition = x;
      xtalcorposition = x-(minCoord+maxCoord)/2;
    }
  // ... and then xtals along the y coordinate.
  else 
    {
      nocorposition = y;
      // And this is slightly more complicated as we have to distinguish between the LAT...
      if ( flightGeom ) 
        {
          minCoord = -1.5*towerPitch + towerPitch*(double)towerIdy - csILength/2;
          maxCoord = -1.5*towerPitch + towerPitch*(double)towerIdy + csILength/2;
        }
      // ...and the CU.
      else 
        {
          minCoord = - csILength/2;
          maxCoord = csILength/2;
        }
      distToEdge = std::min( fabs(y - minCoord), fabs(y - maxCoord) );

      y = fitCentroid.y() + (z - fitCentroid.z()) / fitAxis.z() * fitAxis.y();
      distToEdge = std::min( fabs(y - minCoord), fabs(y - maxCoord) );

      if ( y < minCoord ) 
        {
          y = minCoord;
        }
      if ( y > maxCoord ) 
        {
          y = maxCoord;
        }
      corposition = y;
      xtalcorposition = y-(minCoord+maxCoord)/2;
    }

  if(!recData->longPosCorrected())
    {
      Point corrPos(x,y,z);
      const Event::CalXtalRecData::CalRangeRecData *rngData = recData->getRangeRecData(0); 
      Event::CalXtalRecData::CalRangeRecData *rngDataR = const_cast<Event::CalXtalRecData::CalRangeRecData *> (rngData);
      rngDataR->setPosition(corrPos);
      recData->setLongPosCorrected();
    }

  return xtalcorposition;
}

double MomentsClusterInfo::SetXtalEnergyFromPosition(Event::CalXtalRecData* recData, double xtalposition)
{
  if(!recData->saturated()) return 0;
  if(recData->bothFacesSaturated()) return 0;

  if(recData->energyCorrected()) return 0;

  int tower  = recData->getTower();
  int layer  = recData->getLayer();
  int column = recData->getColumn();

  CalUtil::XtalIdx xtalId(recData->getPackedId());

  CalUtil::AsymType myasymtype = CalUtil::ASYM_SS;
  if(recData->posFaceSaturated() && recData->getRange(0,idents::CalXtalId::NEG)<2) myasymtype = CalUtil::ASYM_LL;
  else if(recData->negFaceSaturated() && recData->getRange(0,idents::CalXtalId::POS)<2) myasymtype = CalUtil::ASYM_LL;

  float myasym;
  m_calCalibSvc->evalAsym(xtalId,myasymtype,(float)xtalposition,myasym);

  float corenergy;
  if(recData->posFaceSaturated())
    {
      corenergy = recData->getEnergy(0,idents::CalXtalId::NEG)*exp(myasym);
      recData->setEnergy(0,idents::CalXtalId::POS,corenergy);
    }
  else
    {
      corenergy = recData->getEnergy(0,idents::CalXtalId::POS)*exp(-myasym);
      recData->setEnergy(0,idents::CalXtalId::NEG,corenergy);
    }
  
  return 0;
}
