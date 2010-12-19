
#include "MomentsClusterInfo.h"
#include "src/Utilities/CalException.h"

#include "CalUtil/CalDefs.h"
#include "Event/Digi/CalDigi.h"

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
MomentsClusterInfo::MomentsClusterInfo(ICalReconSvc* calReconSvc, double transScaleFactor) 
                                     : m_calReconSvc(calReconSvc), 
                                       m_p0(0.,0.,0.),
                                       m_p1(0.,0.,0.),
                                       m_transScaleFactor(transScaleFactor)
{
    m_calnLayers = m_calReconSvc->getCalNLayers();

    m_saturationadc = 4060;

    m_dataVec.clear();

    // Minuit object
    m_minuit = new TMinuit(5);
    
    //Sets the function to be minimized
    m_minuit->SetFCN(fcncal);

    return;
}

/// This makes CalClusters out of associated CalXtalRecData pointers
Event::CalCluster* MomentsClusterInfo::fillClusterInfo(const XtalDataList* xTalVec)
{
    // Create an output cluster
    Event::CalCluster* cluster = 0;

    if (!xTalVec->empty())
    {
        cluster = new Event::CalCluster(0);
        cluster->setProducerName("MomentsClusterInfo") ;

        cluster->clear();

        //  DetectSaturation(calDigiCol);
        DetectSaturation();

        //  getCentroidTransverseInfoOnly(xTalVec);
        fitDirectionCentroid(xTalVec);

        double energy = fillLayerData(xTalVec, cluster);

        fillMomentsData(xTalVec, cluster, energy);
    }

    return cluster;
}

double MomentsClusterInfo::fillLayerData(const XtalDataList* xTalVec, Event::CalCluster* cluster)
{
    //Initialize local variables
    double              ene = 0;                      // Total energy in this cluster
    int                 num_TruncXtals = 0;           // Number of Xtals in the cluster with > 1% of the energy
    Vector              pCluster(0.,0.,0.);           // Cluster position
    std::vector<double> eneLayer(m_calnLayers,0.);    // Energy by layer
    std::vector<Vector> pLayer(m_calnLayers);         // Position by layer
    std::vector<Vector> rmsLayer(m_calnLayers);       // rms by layer

    // Use the fit centroid as the start point
    m_p0 = Point(0.,0.,0.);

    // If the linear fit was performed then use the centroid from that for a first guess
    //if (m_fit_nlayers > 3 && m_fit_zdirection > 0.) m_p0 = Point(m_fit_xcentroid, m_fit_ycentroid, m_fit_zcentroid);

    // Keep track of any energy in crystals with a bad position measurement (resulting from a bad/dead readout)
    double badPosXtalEne = 0.;

    // Compute barycenter and various moments
    // loop over all crystals in the current cluster
    XtalDataList::const_iterator xTalIter;
    for(xTalIter = xTalVec->begin(); xTalIter != xTalVec->end(); xTalIter++)
    {
        // get pointer to the reconstructed data for given crystal
        Event::CalXtalRecData* recData = *xTalIter;
        
        // get reconstructed values
        double eneXtal = recData->getEnergy();                 // crystal energy
        Vector pXtal   = recData->getPosition();               // Vector of crystal position
        int    layer   = (recData->getPackedId()).getLayer();  // layer number

        // Check if bad position measurement
        if (!recData->isPositionGood())
        {
            badPosXtalEne += eneXtal;
            continue;
        }

        // Offset the position to be closer to actual cluster center
        pXtal -= m_p0;

        // update energy of corresponding layer
        eneLayer[layer] += eneXtal;
        
        // update average position of corresponding layer
        Vector ptmp = eneXtal*pXtal;
        pLayer[layer] += ptmp;
        
        // Vector containing squared coordinates, weighted by crystal energy 
        Vector ptmp_sqr(ptmp.x()*pXtal.x(), ptmp.y()*pXtal.y(), ptmp.z()*pXtal.z());
 
        // update quadratic spread, which is proportional to eneXtal;
        // this means, that position error in one crystal
        // is assumed to be 1/sqrt(eneXtal) 
        rmsLayer[layer] += ptmp_sqr;
        
        // update energy sum
        ene  += eneXtal;

        // update cluster position
        pCluster += ptmp;
    }

    // Now take the means
    // if energy sum is not zero - normalize cluster position
    if(ene > 0.) 
    {
        pCluster /= ene;
        pCluster += m_p0;
        cluster->setStatusBit(Event::CalCluster::CENTROID);
    }
    // if energy is zero - set cluster position to non-physical value
    else pCluster = Vector(-1000., -1000., -1000.);
    
    // loop over calorimeter layers
    for(int i = 0; i < m_calnLayers; i++)
    {
        // if energy in the layer is not zero - finalize calculations
        if(eneLayer[i]>0)
        {
            // normalize position in the layer
            pLayer[i] *= (1./eneLayer[i]); 
            
            // normalize quadratic spread in the laye
            rmsLayer[i] *= (1./eneLayer[i]);
            
            // Vector containing the squared average position in each component
            Vector sqrLayer(pLayer[i].x()*pLayer[i].x(),
                            pLayer[i].y()*pLayer[i].y(),
                            pLayer[i].z()*pLayer[i].z());
            
            // the precision of transverse coordinate measurement
            // if there is no fluctuations: 1/sqrt(12) of crystal width
            Vector d;
            double csIWidth = m_calReconSvc->getCalCsIWidth() ;
            if(i%2 == 1) d = Vector(csIWidth*csIWidth/12.,0.,0.);
            else d = Vector(0.,csIWidth*csIWidth/12.,0.);
            
            // subtracting the  squared average position and adding
            // the square of crystal width, divided by 12
            rmsLayer[i] += d-sqrLayer;

            // Reset pLayer to detector coordinates
            pLayer[i] += m_p0;
        }
            
        
        // if energy in the layer is zero - reset position and spread Vectors
        else 
        {
            pLayer[i]   = m_p0;
            rmsLayer[i] = m_p0;
        }


        // Fill the cluster layer data
        Point layerPos(pLayer[i].x(), pLayer[i].y(), pLayer[i].z());

        // Set the data for the vector
        Event::CalClusterLayerData layerData(eneLayer[i], layerPos, rmsLayer[i]);

        cluster->push_back(layerData);
    }

    // Total energy in cluster includes energy in bad position crystals
    ene += badPosXtalEne;

    // loop over all crystals in the current cluster again to compute the number of Truncated Xtals
    for(xTalIter = xTalVec->begin(); xTalIter != xTalVec->end(); xTalIter++)
    {
        // get pointer to the reconstructed data for given crystal
        Event::CalXtalRecData* recData = *xTalIter;
        
        // get reconstructed values
        double eneXtal = recData->getEnergy();                  // crystal energy

        if(eneXtal > .01*ene) num_TruncXtals++;   //NOTE: 1% should be a declared property!!!!!  
    }

    // Correct for Zero Supression using Truncated Xtal count
    ene += .2*num_TruncXtals; //Note .2 MeV/Xtal should be a declared property!!!!!!!!!!!!!

    // Also the Number of Truncated Xtals should be a data member in CalCluster!!!!!!!!!!!

    // Set energy centroid
    Event::CalMomParams momParams(ene, 10*ene, pCluster.x(), pCluster.y(), pCluster.z(),
				  1.,0.,0.,1.,0.,1., 0., 0., 1., 1.,0.,0.,1.,0.,1.);

    // Initial fit parameters
    Event::CalFitParams fitParams(m_fit_nlayers, 0., pCluster.x(), pCluster.y(), pCluster.z(),
				  1.,0.,0.,1.,0.,1., 0., 0., 1., 1.,0.,0.,1.,0.,1.);

    // Use the fit centroid/direction to see moments analysis
    if (m_fit_nlayers > 8 && m_fit_zdirection > 0.)
    {
        fitParams = Event::CalFitParams(m_fit_nlayers, m_fit_chisq, 
                                        m_fit_xcentroid,  m_fit_ycentroid,  m_fit_zcentroid, 
					1.,0.,0.,1.,0.,1.,
                                        m_fit_xdirection, m_fit_ydirection, m_fit_zdirection,
					1.,0.,0.,1.,0.,1.);
    }
    
    // initialize empty CalMSTreeParams container
    Event::CalMSTreeParams treeParams(0.,0.,0,0.,0.,0.,0.,0.,0.);
    // initialize empty prob map - m_classesProb
    std::map <std::string, double> probMap;
    probMap["gam"]=-1;
    
    cluster->initialize(treeParams, fitParams, momParams, probMap, m_Nsaturated, num_TruncXtals);

    return ene;
}

void MomentsClusterInfo::fillMomentsData(const XtalDataList* xTalVec, Event::CalCluster* cluster, double energy)
{
    // Try new utility class
    // Begin by building a Moments Data vector
    m_dataVec.clear();

    Point  centroid = cluster->getMomParams().getCentroid();
    Vector axis     = cluster->getMomParams().getAxis();

    XtalDataList::const_iterator xTalMax = xTalVec->end();

    double rmsDist   = 0.;
    double weightSum = 0.;

    // Loop through the xtals setting the hits to analyze
    for(XtalDataList::const_iterator xTalIter = xTalVec->begin(); xTalIter != xTalVec->end(); xTalIter++)
    {
        if (xTalIter == xTalMax) continue;

        Event::CalXtalRecData* recData = *xTalIter;

        double xPos = recData->getPosition().x();
        double yPos = recData->getPosition().y();
        double zPos = recData->getPosition().z();
    
        // If there are saturated crystals, change the longitudinal position
        if(m_Nsaturated>0 && m_fit_nlayers >= 4 && m_fit_zdirection > 0.)
        {
            if(m_saturated[(int)(recData->getPackedId()).getTower()][(int)(recData->getPackedId()).getLayer()][(int)(recData->getPackedId()).getColumn()])
            {
                Point CorPos = GetCorrectedPosition(recData->getPosition(),
                                    (recData->getPackedId()).getTower(),
                                    (recData->getPackedId()).getLayer(),
                                    (recData->getPackedId()).getColumn());
                xPos = CorPos.x();
                yPos = CorPos.y();
                zPos = CorPos.z();
            }
        }

        Point CrystalPos(xPos, yPos, zPos);

        CalMomentsData momentsData(CrystalPos, recData->getEnergy(), 0.);

        // DC: unused !
        // TU: Actually... bad design on my part but this method must be called. 
        double distToAxis = momentsData.calcDistToAxis(centroid, axis);

        rmsDist   += momentsData.getWeight() * distToAxis * distToAxis;
        weightSum += momentsData.getWeight();
                
        m_dataVec.push_back(momentsData);
    }

    // Get the quick rms of crystals about fit axis
    rmsDist = sqrt(rmsDist / weightSum);

    // Use this to try to remove "isolated" crystals that might otherwise pull the axis
    // Start by sorting the data vector according to distance to axis
    std::sort(m_dataVec.begin(), m_dataVec.end());

    // Look for the point where there is a significant gap in distance to the next xtal
    CalMomentsDataVec::iterator momDataIter;
    double                      envelope = std::min(5.*rmsDist, m_calReconSvc->getCaltowerPitch());
    double                      lastDist = 0.;

    int checkit     = m_dataVec.size();
    int tempCounter = 0;

    // Find the point in the vector where the distance to the axis indicates an outlier
    for(momDataIter = m_dataVec.begin(); momDataIter != m_dataVec.end(); momDataIter++)
    {
        double dist     = (*momDataIter).getDistToAxis();
        double distDiff = dist - envelope;
        double gapDist  = dist - lastDist;

        if (distDiff > 0. && gapDist > 3. * rmsDist) break;

        lastDist = dist;
        tempCounter++;
    }

    // Now remove the outlier crystals
    if (momDataIter != m_dataVec.end()) m_dataVec.erase(momDataIter, m_dataVec.end());

    CalMomentsAnalysis momentsAnalysis;

    double chiSq = momentsAnalysis.doIterativeMomentsAnalysis(m_dataVec, centroid, m_transScaleFactor);

    if (chiSq >= 0.)
    {
        // Get info on iterations
        int nIterations = momentsAnalysis.getNumIterations();
        int nDropped    = momentsAnalysis.getNumDroppedPoints();

        // Get the iterative moments centroid and axis from iterations
        centroid = momentsAnalysis.getMomentsCentroid();
        axis     = momentsAnalysis.getMomentsAxis();

        // Recalculate the moments going back to using all the data points but with
        // the iterated moments centroid
        if (nIterations > 1) chiSq = momentsAnalysis.doMomentsAnalysis(m_dataVec, centroid);
    
        // Extract the values for the moments with all hits present
        double rms_long  = momentsAnalysis.getLongitudinalRms();
        double rms_trans = momentsAnalysis.getTransverseRms();
        double long_asym = momentsAnalysis.getLongAsymmetry();
	double long_skew = momentsAnalysis.getLongSkewness();

        if (!isFinite(rms_long))
        {
            throw CalException("CalMomentsAnalysis computed infinite value for rms_long") ;
        }
        if (!isFinite(rms_trans))
        {
            throw CalException("CalMomentsAnalysis computed infinite value for rms_trans") ;
        }
        if (!isFinite(long_asym))
        {
            throw CalException("CalMomentsAnalysis computed infinite value for long_asym") ;
        }
	if (!isFinite(long_skew))
	{
	    throw CalException("CalMomentsAnalysis computed infinite value for long_skew") ;
        }
    
        int num_TruncXtals = cluster->getNumTruncXtals(); 

        // Store all this information away in the cluster
        //Event::CalParams params(energy, 10*energy,
	//       centroid.x(), centroid.y(), centroid.z(), 1., 0., 0., 1., 0., 1.,
	//       axis.x(),     axis.y(),     axis.z(),     1., 0., 0., 1., 0., 1.);

	// Code for the new CalMomParams class.
	//std::cout << "**********************************************************" << std::endl;
	//std::cout << "CalParams: \n" << params << "\n" << std::endl;
	CLHEP::HepMatrix I_3_3(3, 3, 1);
	Event::CalMomParams momParams (energy, 10*energy, centroid, I_3_3, axis, I_3_3,
				       nIterations, rms_trans, rms_long, long_asym,
				       long_skew, -1.0);
	//std::cout << "CalMomParams: \n" << momParams << std::endl;
	//std::cout << "**********************************************************" << std::endl;

        Event::CalFitParams fitParams(m_fit_nlayers, m_fit_chisq,
                                      m_fit_xcentroid,  m_fit_ycentroid,  m_fit_zcentroid,
				      1., 0., 0., 1., 0., 1.,
                                      m_fit_xdirection, m_fit_ydirection, m_fit_zdirection,
				      1., 0., 0., 1., 0., 1.);

        // initialize empty CalMSTreeParams container - CalMSTreePar
        Event::CalMSTreeParams treeParams(0.,0.,0,0.,0.,0.,0.,0.,0.);
        // initialize empty prob map - m_classesProb
        std::map <std::string, double> probMap;
        probMap["gam"]=-1;

	cluster->initialize(treeParams, fitParams, momParams, probMap, m_Nsaturated,
			    num_TruncXtals);
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
        double eneXtal = recData->getEnergy();                  // crystal energy
        int layer   = (recData->getPackedId()).getLayer();   // layer number
        Vector pXtal   = recData->getPosition();         // Vector of crystal position
        if(layer%2==0)
        {
            myenergy[1] += eneXtal;
            mycentroid[1] += eneXtal*pXtal.y();
        }
        else
        {
            myenergy[0] += eneXtal;
            mycentroid[0] += eneXtal*pXtal.x();
        }
        myenergy[2] += eneXtal;
        mycentroid[2] += eneXtal*pXtal.z();
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
        double eneXtal = recData->getEnergy();                  // crystal energy
        int layer   = (recData->getPackedId()).getLayer();   // layer number
        Vector pXtal   = recData->getPosition();         // Vector of crystal position
        //
        totenergy += eneXtal;
        m_fit_e_layer[layer] += eneXtal;
        mypos[layer][0] += eneXtal*pXtal.x();
        mypos[layer][1] += eneXtal*pXtal.y();
        mypos[layer][2] += eneXtal*pXtal.z();
        mypos2[layer][0] += eneXtal*pXtal.x()*pXtal.x();
        mypos2[layer][1] += eneXtal*pXtal.y()*pXtal.y();
        mypos2[layer][2] += eneXtal*pXtal.z()*pXtal.z();
        //
        if(layer%2==0)
        {
            myenergy[1] += eneXtal;
            mycentroid[1] += eneXtal*pXtal.y();
        }
        else
        {
            myenergy[0] += eneXtal;
            mycentroid[0] += eneXtal*pXtal.x();
        }
        myenergy[2] += eneXtal;
        mycentroid[2] += eneXtal*pXtal.z();
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

Point MomentsClusterInfo::GetCorrectedPosition(Point pcrystal, int itower, int ilayer, int icolumn)
{
    // set the longitudinal position of crystal to the longitudinal position of the extrapolation of the cal direction using only the transverse information
    double xyz[3];
    xyz[0] = pcrystal.x();
    xyz[1] = pcrystal.y();
    xyz[2] = pcrystal.z();
    int itowy = itower/4;
    int itowx = itower-4*itowy;
    double minpos,maxpos;

    if(ilayer%2==0)
    {
        xyz[0] = m_fit_xcentroid+(xyz[2]-m_fit_zcentroid)/m_fit_zdirection*m_fit_xdirection;
        minpos = -1.5*m_calReconSvc->getCaltowerPitch()+m_calReconSvc->getCaltowerPitch()*(double)itowx - m_calReconSvc->getCalCsILength()/2;
        maxpos = -1.5*m_calReconSvc->getCaltowerPitch()+m_calReconSvc->getCaltowerPitch()*(double)itowx + m_calReconSvc->getCalCsILength()/2;
        if(xyz[0]<minpos) xyz[0] = minpos;
        if(xyz[0]>maxpos) xyz[0] = maxpos;
    }
    else
    {
        xyz[1] = m_fit_ycentroid+(xyz[2]-m_fit_zcentroid)/m_fit_zdirection*m_fit_ydirection;
        if(m_calReconSvc->getCalFlightGeom())
        {
            minpos = -1.5*m_calReconSvc->getCaltowerPitch()+m_calReconSvc->getCaltowerPitch()*(double)itowy - m_calReconSvc->getCalCsILength()/2;
            maxpos = -1.5*m_calReconSvc->getCaltowerPitch()+m_calReconSvc->getCaltowerPitch()*(double)itowy + m_calReconSvc->getCalCsILength()/2;
        }
        else
        {
            minpos = - m_calReconSvc->getCalCsILength()/2;
            maxpos = m_calReconSvc->getCalCsILength()/2;
        }
        if(xyz[1]<minpos) xyz[1] = minpos;
        if(xyz[1]>maxpos) xyz[1] = maxpos;
    }

    Point CorPos(xyz[0],xyz[1],xyz[2]);

    return CorPos;
}

int MomentsClusterInfo::DetectSaturation()
{
    m_Nsaturated = 0;
    std::memset(m_saturated, false, 16*8*12*sizeof(bool));

    const Event::CalDigiCol *calDigiCol = m_calReconSvc->getDigis();

    if(calDigiCol==NULL) return 0;

    for (Event::CalDigiCol::const_iterator digiIter = calDigiCol->begin(); digiIter != calDigiCol->end(); digiIter++)
    {
        const Event::CalDigi calDigi = **digiIter;
        idents::CalXtalId    id      = calDigi.getPackedId();

        CalUtil::XtalIdx xtalIdx(calDigi.getPackedId());

        for (Event::CalDigi::CalXtalReadoutCol::const_iterator ro =  calDigi.getReadoutCol().begin();ro != calDigi.getReadoutCol().end();ro++)
        {
            float adcP     = ro->getAdc(idents::CalXtalId::POS);
            float adcN     = ro->getAdc(idents::CalXtalId::NEG);
            int   rangepos = ro->getRange(idents::CalXtalId::POS);
            int   rangeneg = ro->getRange(idents::CalXtalId::NEG);

            if( (rangepos==3 && adcP>=m_saturationadc) || (rangeneg==3 && adcN>=m_saturationadc))
            {
                m_saturated[id.getTower()][id.getLayer()][id.getColumn()] = true;
                ++m_Nsaturated;
            }
        }
    }
  
    return 0;
}
