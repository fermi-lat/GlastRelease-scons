
#include "MomentsClusterInfo.h"
#include "src/Utilities/CalException.h"

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
 
/// Constructor
MomentsClusterInfo::MomentsClusterInfo(const ICalReconSvc* calReconSvc, double transScaleFactor) 
                                     : m_calReconSvc(calReconSvc), 
                                       m_p0(0.,0.,0.),
                                       m_transScaleFactor(transScaleFactor)
{
    m_calnLayers = m_calReconSvc->getCalNLayers();

    m_dataVec.clear();

    return;
}

/// This makes CalClusters out of associated CalXtalRecData pointers
Event::CalCluster* MomentsClusterInfo::fillClusterInfo(const XtalDataVec* xTalVec)
{
    // Create an output cluster
    Event::CalCluster* cluster = 0;
    
    if (!xTalVec->empty())
    {
        cluster = new Event::CalCluster(0);
        cluster->setProducerName("MomentsClusterInfo") ;

        cluster->clear();

        double energy = fillLayerData(xTalVec, cluster);

        fillMomentsData(xTalVec, cluster, energy);
    }

    return cluster;
}

double MomentsClusterInfo::fillLayerData(const XtalDataVec* xTalVec, Event::CalCluster* cluster)
{
    //Initialize local variables
    double ene = 0;                                   // Total energy in this cluster
    int num_TruncXtals = 0;                           // Number of Xtals in the cluster with > 1% of the energy
    Vector pCluster(0.,0.,0.);                        // Cluster position
    std::vector<double> eneLayer(m_calnLayers,0.);    // Energy by layer
    std::vector<Vector> pLayer(m_calnLayers);         // Position by layer
    std::vector<Vector> rmsLayer(m_calnLayers);       // rms by layer

    // Compute barycenter and various moments
    // loop over all crystals in the current cluster
    XtalDataVec::const_iterator xTalIter;
    for(xTalIter = xTalVec->begin(); xTalIter != xTalVec->end(); xTalIter++)
    {
        // get pointer to the reconstructed data for given crystal
        Event::CalXtalRecData* recData = *xTalIter;
        
        // get reconstructed values
        double eneXtal = recData->getEnergy();                  // crystal energy
        Vector pXtal   = recData->getPosition() - m_p0;         // Vector of crystal position
        int    layer   = (recData->getPackedId()).getLayer();   // layer number

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
        cluster->setStatusBit(Event::CalCluster::CENTROID);
    }
    // if energy is zero - set cluster position to non-physical value
    else pCluster = Vector(-1000., -1000., -1000.);
    
    // loop over calorimeter layers
    int i ;
    for( i = 0; i < m_calnLayers; i++)
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
    Event::CalParams params(ene, 10*ene, pCluster.x(), pCluster.y(), pCluster.z(), 1.,0.,0.,1.,0.,1.,
                                               0., 0., 1.,   1.,0.,0.,1.,0.,1.);
    cluster->initialize(params, 0., 0., 0., num_TruncXtals);

    return ene;
}

void MomentsClusterInfo::fillMomentsData(const XtalDataVec* xTalVec, Event::CalCluster* cluster, double energy)
{
    // Try new utility class
    // Begin by building a Moments Data vector
    m_dataVec.clear();

    Point  centroid = cluster->getCalParams().getCentroid();
    Vector axis     = cluster->getCalParams().getAxis();

    XtalDataVec::const_iterator xTalMax = xTalVec->end();

    // Loop through the xtals setting the hits to analyze
    for(XtalDataVec::const_iterator xTalIter = xTalVec->begin(); xTalIter != xTalVec->end(); xTalIter++)
    {
        if (xTalIter == xTalMax) continue;

        Event::CalXtalRecData* recData = *xTalIter;

        CalMomentsData momentsData(recData->getPosition(), recData->getEnergy(), 0.);

        // DC: unused !
        //double distToAxis = momentsData.calcDistToAxis(centroid, axis);
                
        m_dataVec.push_back(momentsData);
    }

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
    
        int num_TruncXtals = cluster->getNumTruncXtals(); 

        // Store all this information away in the cluster
        Event::CalParams params(energy, 10*energy, centroid.x(), centroid.y(), centroid.z(), 1.,0.,0.,1.,0.,1.,
                                                   axis.x(),     axis.y(),     axis.z(),     1.,0.,0.,1.,0.,1.);

        cluster->initialize(params, rms_long, rms_trans, long_asym, num_TruncXtals);
        cluster->setStatusBit(Event::CalCluster::MOMENTS); 
    }

    return;
}
