/**
 * @class SimpleClusterTool
 *
 * @brief Implements a Gaudi Tool for performing very simple clustering in the Cal 
 *
 * @author Tracy Usher
 *
 * File and Version Information:
 *      $Header$
 */

// Tool and Gaudi related stuff
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/DeclareFactoryEntries.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "ICluster.h"

// TDS related stuff
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"

typedef  std::vector<Event::CalXtalRecData*> xTalDataVec;

class SimpleClusterTool : public AlgTool, virtual public ICluster
{
public:
    /// Standard Gaudi Tool interface constructor
    SimpleClusterTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~SimpleClusterTool();

	/// @brief Intialization of the tool
    StatusCode initialize();

    /// @brief Execute is defined by ICluster
    StatusCode execute();

    /// @brief Finalization of the tool
    StatusCode finalize();

    /// @brief This does the cluster finding
    StatusCode findClusters(Event::CalXtalRecCol* calXtalRecCol);

    /// @brief This initializes the data set
    void setClusterCol(Event::CalClusterCol* calClusterCol) {m_calClusterCol = calClusterCol;};

protected:
    Vector     Fit_Direction(std::vector<Vector> pos, std::vector<Vector> sigma2, int nlayers);

private:
    /// This finds the next highest energy cluster in a vector of CalXtalRecData pointers
    xTalDataVec            getCluster(xTalDataVec& xTalVec);
    xTalDataVec            getXtalsInLayer(xTalDataVec& xTalVec, Event::CalXtalRecData* xTal);
    Event::CalXtalRecData* getNearestXtalInSameLayer(xTalDataVec& xTalVec, xTalDataVec& NNvec, Event::CalXtalRecData* xTal);
    Event::CalXtalRecData* getNearestXtalInDiffLayer(xTalDataVec& xTalVec, Event::CalXtalRecData* xTal, int layer);

    /// This makes a CalCluster out of associated CalXtalRecData pointers
    void               saveCluster(xTalDataVec& xTalVec);

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc*    m_dataSvc;

    /// pointer to GlasDetSvc
    IGlastDetSvc* detSvc;
       
	//! crystal width
    double m_CsIWidth;

    //! crystal height
    double m_CsIHeight;

    //! number of layers
    int m_CalnLayers;

	//! reconstructed data for crystals, the input of the reconstruction
	Event::CalXtalRecCol* m_calXtalRecCol;
		
	Event::CalClusterCol* m_calClusterCol;
};

//static ToolFactory<SimpleClusterTool> s_factory;
//const IToolFactory& SimpleClusterToolFactory = s_factory;
DECLARE_TOOL_FACTORY(SimpleClusterTool) ;

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

SimpleClusterTool::SimpleClusterTool(const std::string& type, const std::string& name, const IInterface* parent) :
                    AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ICluster>(this);

    return;
}

// 
// Cleanup memory on exit
//
SimpleClusterTool::~SimpleClusterTool()
{
    return;
}
//
// Initialization of the tool here
//

StatusCode SimpleClusterTool::initialize()
{	
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
	log << MSG::INFO << "BEGIN initialize()" <<endreq;

    //Set the properties
    setProperties();

    //Locate and store a pointer to the data service
    if( (sc = service("EventDataSvc", m_dataSvc)).isFailure() ) 
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    IParticlePropertySvc*  partPropSvc;
    if( (sc = service("ParticlePropertySvc", partPropSvc,true)).isFailure() ) 
    {
        throw GaudiException("Service [ParticlePropertySvc] not found", name(), sc);
    }
    // get pointer to GlastDetSvc
    sc = service("GlastDetSvc", detSvc);
    
    // if GlastDetSvc isn't available - put error message and return
    if(sc.isFailure())
    {
        log << MSG::ERROR << "GlastDetSvc could not be found" <<endreq;
        return sc;
    }
    
    // extracting detector geometry constants from xml file
    double value;
    if(!detSvc->getNumericConstByName(std::string("CALnLayer"), &value)) 
    {
        log << MSG::ERROR << " constant " << " CALnLayer "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } else m_CalnLayers = int(value);
    
    if(!detSvc->getNumericConstByName(std::string("CsIWidth"),&m_CsIWidth))
    {
        log << MSG::ERROR << " constant " << " CsIWidth "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } 
    
    if(!detSvc->getNumericConstByName(std::string("CsIHeight"),&m_CsIHeight))
    {
        log << MSG::ERROR << " constant " << " CsIHeight "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } 

	log << MSG::INFO << "END initialize()" <<endreq;
    return sc;
}


StatusCode SimpleClusterTool::finalize()
{
	StatusCode sc = StatusCode::SUCCESS;

	return sc;
}

StatusCode SimpleClusterTool::execute()
{
	StatusCode sc = StatusCode::SUCCESS;

	return sc;
}


Vector SimpleClusterTool::Fit_Direction(std::vector<Vector> pos,
                                        std::vector<Vector> sigma2,
                                        int nlayers)
//
// Purpose and Method:
//       find the particle direction from average positions in each
//       layer and their quadratic errors. The fit is made independently
//       in XZ and YZ planes for odd and even layers, respectively.
//       Then the 3-vector of particle direction is calculated from these 
//       2 projections.
//       Only position information based on the transversal crystal position
//       is used. The position along the crystal, calculated from signal
//       asymmetry is not used.
//
// Inputs:
//         pos      - average position for each calorimeter layer
//         sigma2   - quadratic error of position measurement for each layer
//                    in each direction  
//         nlayers  - number of calorimeter layers
//
// Returned value:    3-Vector of reconstructred particle direction
//
                                     
{
    
    MsgStream log(msgSvc(), name());
    
    // sigma2.z() is useless here no matter its value.
    double cov_xz = 0;  // covariance x,z
    double cov_yz = 0;  // covariance y,z
    double mx=0;        // mean x
    double my=0;        // mean y
    double mz1=0;       // mean z for x pos
    double mz2=0;       // mean z for y pos
    double norm1=0;     // sum of weights for odd layers
    double norm2=0;     // sum of weights for even layers
    double var_z1=0;    // variance of z for odd layers	
    double var_z2=0;    // variance of z for even layers
    
    // number of layers with non-zero energy deposition
    // in X and Y direction, respectively
    int nlx=0,nly=0;
    
    
    // "non-physical vector of direction, which is returned
    // if fit is imposible due to insufficient number of hitted layers
    Vector nodir(-1000.,-1000.,-1000.);
    
    
    // loop over calorimeter layers
    for(int il=0;il<m_CalnLayers;il++)
    {                
        // For the moment forget about longitudinal position
        
        // odd layers used for XZ fit
        if(il%2==1)
        {
            
            // only if the spread in X direction is not zero,
            // which is the case if there is non-zero energy
            // deposition in this layer
            if (sigma2[il].x()>0.)
            {
                nlx++; // counting layers used for the fit in XZ plane 
                
                // calculate weighting coefficient for this layer
                double err = 1/sigma2[il].x(); 
                
                // calculate sums for least square linear fit in XZ plane
                cov_xz += pos[il].x()*pos[il].z()*err;
                var_z1 += pos[il].z()*pos[il].z()*err;
                mx += pos[il].x()*err;
                mz1 += pos[il].z()*err;
                norm1 += err;
            }
        }
        // even layers used for YZ fit
        else
        {
            // only if the spread in Y direction is not zero,
            // which is the case if there is non-zero energy
            // deposition in this layer
            if(sigma2[il].y()>0.)
            {
                
                nly++; // counting layers used for the fit in YZ plane 
                
                // calculate weighting coefficient for this layer
                double err = 1/sigma2[il].y();
                
                
                // calculate sums for least square linear fit in YZ plane
                cov_yz += pos[il].y()*pos[il].z()*err;
                var_z2 += pos[il].z()*pos[il].z()*err;
                my += pos[il].y()*err;
                mz2 += pos[il].z()*err;
                norm2 += err;
            }
        }
    }		
    
    // linear fit requires at least 2 hitted layers in both XZ and YZ planes
    // otherwise non-physical direction is returned
    // which means that direction hasn't been found
    if(nlx <2 || nly < 2 )return nodir;
    
    
    
    
    mx /= norm1;
    mz1 /= norm1;
    cov_xz /= norm1;
    cov_xz -= mx*mz1;
    var_z1 /= norm1;
    var_z1 -= mz1*mz1;
    
    // protection against dividing by 0 in the next statment
    if(var_z1 == 0) return nodir;
    
    // Now we have cov(x,z) and var(z) we can
    // deduce slope in XZ plane
    double tgthx = cov_xz/var_z1;
    
    
    my /= norm2;
    mz2 /= norm2;
    cov_yz /= norm2;
    cov_yz -= my*mz2;
    var_z2 /= norm2;
    var_z2 -= mz2*mz2;
    
    // protection against dividing by 0 in the next statment
    if(var_z2 == 0) return nodir;
    
    // Now we have cov(y,z) and var(z) we can
    // deduce slope in YZ plane
    double tgthy = cov_yz/var_z2;
    
    // combining slope in XZ and YZ planes to get normalized 3-vector
    // of particle direction
    double tgtheta_sqr = tgthx*tgthx+tgthy*tgthy;
    double costheta = 1/sqrt(1+tgtheta_sqr);
    Vector dir(costheta*tgthx,costheta*tgthy,costheta);
    return dir;
}


StatusCode SimpleClusterTool::findClusters(Event::CalXtalRecCol* calXtalRecCol)

//Purpose and method:
//
//   This function performs the calorimeter cluster reconstruction.
//   The main actions are:
//      - calculate energy sum
//                  energy per layer
//                  average position per layer
//                  quadratic spread per layer
//      - fit the particle direction using Fit_Direction() function
//      - store all calculated quantities in CalCluster objects
// 
// TDS input: CalXtalRecCol
// TDS output: CalClustersCol


{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    //Copy pointers to crystal objects for local use
    xTalDataVec xTalData;
    xTalData.clear();

    Event::CalXtalRecCol::const_iterator it;
    for (it = calXtalRecCol->begin(); it != calXtalRecCol->end(); it++)
    {
        // get pointer to the reconstructed data for given crystal
		Event::CalXtalRecData* recData = *it;
        xTalData.push_back(recData);
    }

    //Make clusters
    if (int numLeft  = xTalData.size() > 0)
    {
        while(numLeft > 0)
        {
            xTalDataVec cluster = getCluster(xTalData);

            int numXtals = cluster.size();

            numLeft  = xTalData.size();

            saveCluster(cluster);
        }
    }
    // Always store a zero cluster so downstream code runs ok
    else
    {
        saveCluster(xTalData);
    }

	return sc;
}

/// This finds the next highest energy cluster in a vector of CalXtalRecData pointers
xTalDataVec SimpleClusterTool::getCluster(xTalDataVec& xTalVec)
{
    xTalDataVec cluster;
    cluster.clear();

    //Start by finding the highest energy crystal in the vector
    double bestE = 0.;
    xTalDataVec::iterator bestIter;
    xTalDataVec::iterator xTalVecIter;
    for(xTalVecIter = xTalVec.begin(); xTalVecIter != xTalVec.end(); xTalVecIter++)
    {
        Event::CalXtalRecData* xTalRec = *xTalVecIter;

        if (xTalRec->getEnergy() > bestE)
        {
            bestIter = xTalVecIter;
            bestE    = xTalRec->getEnergy();
        }
    }

    //Add this crystal to our new cluster list and remove from the old list
    Event::CalXtalRecData* xTal = *bestIter;
    cluster.push_back(*bestIter);
    xTalVec.erase(bestIter);

    //Now build up a list of all xTals in this layer which are neighbors 
    Event::CalXtalRecData* nxtXtal = getNearestXtalInSameLayer(xTalVec, cluster, xTal);

    if (nxtXtal) cluster.push_back(nxtXtal);

    //Do again to pick up neighbor on the other side (we need a better algorithm here...) 
    nxtXtal = getNearestXtalInSameLayer(xTalVec, cluster, xTal);

    if (nxtXtal) cluster.push_back(nxtXtal);

    //Get the current xTal id
    const idents::CalXtalId xTalId = xTal->getPackedId();

    //Loop "up" (if possible) associating crystals above us
    Event::CalXtalRecData* bestXtal = xTal;
    for(int layer = xTalId.getLayer() - 1; layer >= 0; layer--)
    {
        //Finds the closest xTal in the next layer, if successful removes it from list
        bestXtal = getNearestXtalInDiffLayer(xTalVec, bestXtal, layer);

        //No crystal signifies we are done
        if (bestXtal == 0) break;
        //if (bestXtal == 0) continue;

        //Add to cluster
        cluster.push_back(bestXtal);

        //Searches for nearest neigbors, removes from current list when found
        nxtXtal = getNearestXtalInSameLayer(xTalVec, cluster, bestXtal);
        if (nxtXtal) cluster.push_back(nxtXtal);
        nxtXtal = getNearestXtalInSameLayer(xTalVec, cluster, bestXtal);
        if (nxtXtal) cluster.push_back(nxtXtal);
    }

    //Loop "down" (if possible) associating xTals below us
    bestXtal = xTal;
    for(int layer = xTalId.getLayer() + 1; layer < m_CalnLayers; layer++)
    {
        bestXtal = getNearestXtalInDiffLayer(xTalVec, bestXtal, layer);

        if (bestXtal == 0) break;
        //if (bestXtal == 0) continue;

        cluster.push_back(bestXtal);

        nxtXtal = getNearestXtalInSameLayer(xTalVec, cluster, bestXtal);
        if (nxtXtal) cluster.push_back(nxtXtal);
        nxtXtal = getNearestXtalInSameLayer(xTalVec, cluster, bestXtal);
        if (nxtXtal) cluster.push_back(nxtXtal);
    }

    // Finished! 
    return cluster;
}

xTalDataVec SimpleClusterTool::getXtalsInLayer(xTalDataVec& xTalVec, Event::CalXtalRecData* xTal)
{
    xTalDataVec newVec;
    newVec.clear();

    //Extract the ID information we need from current crystal
    const idents::CalXtalId xTalId   = xTal->getPackedId();
    int                     curTower = xTalId.getTower();
    int                     curLayer = xTalId.getLayer();

    //Loop through input vector of Xtals looking for a match to this layer
    xTalDataVec::iterator xTalVecIter;
    for(xTalVecIter = xTalVec.begin(); xTalVecIter != xTalVec.end(); xTalVecIter++)
    {
        Event::CalXtalRecData* xTalRec = *xTalVecIter;

        const idents::CalXtalId xTalId  = xTalRec->getPackedId();

        //If same tower and layer then store away
        if (xTalId.getTower() == curTower && xTalId.getLayer() == curLayer) newVec.push_back(xTalRec);
    }

    return newVec;
}

Event::CalXtalRecData* SimpleClusterTool::getNearestXtalInSameLayer(xTalDataVec& xTalVec, xTalDataVec& NNvec, Event::CalXtalRecData* xTal)
{
    //Extract the ID information we need from current crystal
    const idents::CalXtalId xTalId    = xTal->getPackedId();
    int                     curTower  = xTalId.getTower();
    int                     curLayer  = xTalId.getLayer();
    int                     curColumn = xTalId.getColumn();

    //Loop through input vector of Xtals looking for a nearest neighbor match
    xTalDataVec::iterator xTalVecIter;
    for(xTalVecIter = xTalVec.begin(); xTalVecIter != xTalVec.end(); xTalVecIter++)
    {
        Event::CalXtalRecData* xTalRec  = *xTalVecIter;
        const idents::CalXtalId xTalId  = xTalRec->getPackedId();
        int                     column  = xTalId.getColumn(); 
        int                     colDiff = curColumn - column;

        //Only accept Xtal if it is exactly next to the current one
        if (curTower == xTalId.getTower() && curLayer == xTalId.getLayer() && abs(colDiff) < 2)
        {
            //Remove this Xtal from the current list
            xTalVec.erase(xTalVecIter);

            //Look for the nearest neighbor to this Xtal
            Event::CalXtalRecData* nextXtal = getNearestXtalInSameLayer(xTalVec, NNvec, xTalRec);

            //If one found, add to the Nearest Neighbor list
            if (nextXtal) NNvec.push_back(nextXtal);

            //return the current Xtal
            return xTalRec;
        }
    }

    //If we got here then nothing found, return a null pointer
    return 0;
}

Event::CalXtalRecData* SimpleClusterTool::getNearestXtalInDiffLayer(xTalDataVec& xTalVec, Event::CalXtalRecData* xTal, int layer)
{
    Event::CalXtalRecData* newXtal = 0;

    //Extract the ID information we need from current crystal, also get position
    const idents::CalXtalId xTalId   = xTal->getPackedId();
    int                     curTower = xTalId.getTower();
    Point                   curPos   = xTal->getPosition();
    bool                    isXlyr   = layer % 2 == 0;

    //Keep track of the best match
    xTalDataVec::iterator bestXtalIter = xTalVec.end();
    double                bestDist     = 10.;

    //Loop through input vector of Xtals looking for a nearest neighbor match
    xTalDataVec::iterator xTalVecIter;
    for(xTalVecIter = xTalVec.begin(); xTalVecIter != xTalVec.end(); xTalVecIter++)
    {
        Event::CalXtalRecData*  xTalRec = *xTalVecIter;
        const idents::CalXtalId xTalId  = xTalRec->getPackedId();

        //Only accept Xtal if it is exactly next to the current one
        if (curTower == xTalId.getTower() && layer == xTalId.getLayer())
        {
            //Compute distance to this xTal
            Vector distVec = curPos - xTalRec->getPosition();
            double dist    = distVec.magnitude() / m_CsIHeight;

            if (dist < bestDist)
            {
                bestDist     = dist;
                bestXtalIter = xTalVecIter;
            }
        }
    }

    //Did we get an acceptable crystal?
    if (bestXtalIter != xTalVec.end())
    {
        newXtal = *bestXtalIter;
        xTalVec.erase(bestXtalIter);
    }

    //Return out find...
    return newXtal;
}

/// This makes a CalCluster out of associated CalXtalRecData pointers
void SimpleClusterTool::saveCluster(xTalDataVec& xTalVec)
{
    const Point p0(0.,0.,0.);  

    //Initialize variables
    double ene = 0;                                 // Total energy in this cluster
    Vector pCluster(0.,0.,0.);                      // Cluster position
    std::vector<double> eneLayer(m_CalnLayers);     // Energy by layer
    std::vector<Vector> pLayer(m_CalnLayers);       // Position by layer
    std::vector<Vector> rmsLayer(m_CalnLayers);     // rms by layer

    //Zero the energy by layer 
    for(int i = 0; i < m_CalnLayers; i++)
    {
        eneLayer[i] = 0.;
    }


    // Compute barycenter and various moments
    
    // loop over all crystals in the current cluster
    xTalDataVec::iterator xTalIter;
    for(xTalIter = xTalVec.begin(); xTalIter != xTalVec.end(); xTalIter++)
    {
        // get pointer to the reconstructed data for given crystal
		Event::CalXtalRecData* recData = *xTalIter;
        
        // get reconstructed values
        double eneXtal = recData->getEnergy();                // crystal energy
        Vector pXtal   = recData->getPosition() - p0;         // Vector of crystal position
        int    layer   = (recData->getPackedId()).getLayer(); // layer number

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
    if(ene > 0.) pCluster /= ene; 
 	// if energy is zero - set cluster position to non-physical value
    else pCluster = Vector(-1000., -1000., -1000.);
    
    // loop over calorimeter layers
    for(int i = 0; i < m_CalnLayers; i++)
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
            if(i%2 == 1) d = Vector(m_CsIWidth*m_CsIWidth/12.,0.,0.);
            else d = Vector(0.,m_CsIWidth*m_CsIWidth/12.,0.);
            
            // subtracting the  squared average position and adding
            // the square of crystal width, divided by 12
            rmsLayer[i] += d-sqrLayer;
        }
            
        
        // if energy in the layer is zero - reset position and spread Vectors
        else 
        {
            pLayer[i]=p0;
            rmsLayer[i]=p0;
        }
    }
    
    // Now sum the different rms to have one transverse and one longitudinal rms
    double rms_trans = 0;
    double rms_long = 0;
    std::vector<Vector> posrel(m_CalnLayers);
    
    for(int ilayer = 0; ilayer < m_CalnLayers; ilayer++)
    {
        posrel[ilayer]=pLayer[ilayer]-pCluster;
       
        // Sum alternatively the rms
        if(ilayer%2)
        {
            rms_trans += rmsLayer[ilayer].y();
            rms_long += rmsLayer[ilayer].x();
        }
        else
        {
            rms_trans += rmsLayer[ilayer].x();
            rms_long += rmsLayer[ilayer].y();
        }
    }
 
    // Compute direction using the positions and rms per layer
    Vector caldir = Fit_Direction(pLayer, rmsLayer, m_CalnLayers);

    // Fill CalCluster data
    Event::CalCluster* cl = new Event::CalCluster(ene, pCluster + p0);

    cl->initialize(ene, eneLayer, pLayer, rmsLayer, rms_long, rms_trans, caldir, 0.);

    m_calClusterCol->add(cl);

    return;
}
