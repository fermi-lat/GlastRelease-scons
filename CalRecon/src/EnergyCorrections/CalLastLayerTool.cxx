#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"

#include <CalRecon/ICalEnergyCorr.h>
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

// to access an XML containing Digi parameters file
#include "xmlBase/IFile.h"

class CalLastLayerTool  : public AlgTool, virtual public ICalEnergyCorr
{
public:
    
    //! destructor
    CalLastLayerTool( const std::string& type, const std::string& name, const IInterface* parent);
     ~CalLastLayerTool() {}; 
    
     StatusCode initialize();

        
//! Leakage correction using last layer correlation
/*!This method uses the correlation between the escaping energy and
*  the energy deposited in the last layer of the calorimeter.  Indeed,
* the last layer carries the most important information concerning the leaking
* energy:  the total number of particles escaping through the back should be
* nearly proportional to the  energy deposited in the last layer.
*  The measured signal in that layer can therefore be modified to account
* for the leaking energy. 
* We used the Monte Carlo simulation of the GLAST beam test configuration to 
* determine this correlation at several energies, from 2 GeV up to 40 GeV.  
* For one particular incident energy, the bidimensionnal distribution of 
* the energy escaping and the energy deposited in the last layer can be 
* fitted by a simple linear function:
* \f[ E_{leak} = \alpha E_{last} \f]
* The coefficient \f$\alpha\f$   is parametrized as a function of logarithm
 of incident energy: 
*\f[ \alpha = (p_0 + p_1 * lnE )/(1+exp(-p_2*(lnE - p_3))) \f]
* where coefficients \f$ p_i \f$ are the linear functions
* of \f$ cos(\theta) \f$:
* \f[ p_i = a_i + b_i * cos(\theta) \f]
* The coefficients \f$ a_i \f$ and \f$ b_i \f$ are obtained by fitting MC
* simulation results. 
* The value of \f$ cos(\theta) \f$ is determined from tracker reconstruction
* or, if it is not available, from fitting calorimeter data
* (see Fit_Direction() function for details).  
* The only information on the incident energy we have initially is       
*  the total measured energy
* \f$E_m\f$, thus we have to use it as the estimator of incident energy \f$E_0\f$.
*  The reconstructed energy is then:
* \f[ E_{rec} = E_m + \alpha(E_m) E_{last} \f]
* To improve the result, we make one iteration using the new energy estimator to determine
* the  correct value of \f$\alpha\f$.  
*
*\par The method takes 2 arguments:
*\param eTotal Total energy measured in the calorimeter in MeV
*\param elast  Total energy measured in the last layer in MeV
*
*\return Corrected energy in MeV
*
* \warning Give biased results when the last layer gains are misaligned
*
* \author Regis Terrier
*
* \b Revision: 
* - 10/17/00    RT    comments added
* - 09/15/00    RT    bias correction added, tuned on MC data
* - 08/20/00    RT    first implementation
*/
           
     // worker function to calculate energy correction 
     Event::CalCorToolResult* doEnergyCorr(Event::CalCluster*, Event::TkrVertex* ) ;

    
 private:

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc* m_dataSvc;

    /// Detector Service
    IGlastDetSvc *    m_detSvc; 

    /// input XML file containing parameters for Digitization
    std::string m_xmlFile;
    
    //! CAL number of layers
    int m_calNLayers ;
    //! CAL crystal width
    double m_calCsIWidth ;
    //! CAL crystal height
    double m_calCsIHeight ;
    
/// correlation factor 1 
    float m_c0;    
/// correlation factor 2
    float m_c1;    
/// correlation factor 3
    float m_c2;    
/// correlation factor 4
    float m_c3;    
    
/// bias factor 1 
    float m_b0;
/// bias factor 2     
    float m_b1;    

/// golp factor 0 
    float m_k0;
/// golp factor 1 
    float m_k1;
/// golp factor 2     
    float m_k2;    

};

#include "GaudiKernel/DeclareFactoryEntries.h"
DECLARE_TOOL_FACTORY(CalLastLayerTool) ;

CalLastLayerTool::CalLastLayerTool( const std::string& type, 
                                      const std::string& name, 
                                      const IInterface* parent)
                                    : AlgTool(type,name,parent)
{

    // declare base interface for all consecutive concrete classes
    declareInterface<ICalEnergyCorr>(this);
    declareProperty ("xmlFile", m_xmlFile="$(CALRECONROOT)/xml/CalLayer.xml");
};


StatusCode CalLastLayerTool::initialize()
{
    // This function does following initialization actions:
    //    - extracts geometry constants from xml file using GlastDetSvc
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    log << MSG::INFO << "Initializing CalLastLayerTool" <<endreq;

    //Locate and store a pointer to the data service which allows access to the TDS
    if ((sc = service("EventDataSvc", m_dataSvc)).isFailure())
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    if ((sc = service("GlastDetSvc", m_detSvc, true)).isFailure())
    { 
        throw GaudiException("Service [GlastDetSvc] not found", name(), sc);
    }
    
    // extracting detector geometry constants from xml file
    if (!m_detSvc->getNumericConstByName(std::string("CALnLayer"), &m_calNLayers)) 
    { 
        throw GaudiException("GlastDetSvc cannot find [CALnLayer]", name(), StatusCode::FAILURE);
    }
  
    if (!m_detSvc->getNumericConstByName(std::string("CsIWidth"),&m_calCsIWidth))
    { 
        throw GaudiException("GlastDetSvc cannot find [CsIWidth]", name(), StatusCode::FAILURE);
    }
  
    if (!m_detSvc->getNumericConstByName(std::string("CsIHeight"),&m_calCsIHeight))
    { 
        throw GaudiException("GlastDetSvc cannot find [CsIHeight]", name(), StatusCode::FAILURE);
    }

 // Read in the parameters from the XML file
    xmlBase::IFile m_ifile(m_xmlFile.c_str());
    if (m_ifile.contains("lastlayer","c0") && m_ifile.contains("lastlayer","c1") && m_ifile.contains("lastlayer","c2") && m_ifile.contains("lastlayer","c3") && m_ifile.contains("lastlayer","b0") && m_ifile.contains("lastlayer","b1") ) 
    {
        m_c0 = m_ifile.getDouble("lastlayer", "c0");
        log << MSG::INFO << " value for c0 " << m_c0 << endreq;
        m_c1 = m_ifile.getDouble("lastlayer", "c1");
        log << MSG::INFO << " value for c1 " << m_c1 << endreq;
        m_c2 = m_ifile.getDouble("lastlayer", "c2");
        log << MSG::INFO << " value for c0 " << m_c2 << endreq;
        m_c3 = m_ifile.getDouble("lastlayer", "c3");
        log << MSG::INFO << " value for c3 " << m_c3 << endreq;
        m_b0 = m_ifile.getDouble("lastlayer", "b0");
        log << MSG::INFO << " value for b0 " << m_b0 << endreq;
        m_b1 = m_ifile.getDouble("lastlayer", "b1");
        log << MSG::INFO << " value for b1 " << m_b1 << endreq;
        m_k0 = m_ifile.getDouble("lastlayer", "k0");
        log << MSG::INFO << " value for k0 " << m_k0 << endreq;
        m_k1 = m_ifile.getDouble("lastlayer", "k1");
        log << MSG::INFO << " value for k1 " << m_k1 << endreq;
        m_k2 = m_ifile.getDouble("lastlayer", "k2");
        log << MSG::INFO << " value for k2 " << m_k2 << endreq;
    }
    else return StatusCode::FAILURE;
     
    return sc;
}


//Purpose and method:
//
//   This function performs 
//   The main actions are:
//      - calculate particle energy by last layer correction method
//          using Leak() function
// 
// TDS input: CalCluster
// TDS output: CalClusters

Event::CalCorToolResult* CalLastLayerTool::doEnergyCorr(Event::CalCluster* cluster, Event::TkrVertex* vertex)
{
    Event::CalCorToolResult * corResult = new Event::CalCorToolResult();
    corResult->setStatusBit(Event::CalCorToolResult::VALIDPARAMS);
    corResult->setCorrectionName(type());
    corResult->setChiSquare(1.);
    
    Event::CalParams params = cluster->getCalParams() ;

    double eTotal = params.getEnergy() ;
    
    double llStatus = 0 ;
  
    if (eTotal<500.) 
    {
        llStatus = -1. ;
    }
    else 
    {
        // if reconstructed tracker data doesn't exist - put the debugging message
        if (vertex == 0) 
        {
            llStatus = -2 ;
            //log << MSG::DEBUG << "No TKR Reconstruction available " << endreq;
            //return sc;
        }
        // if data exist and number of tracks not zero 
        // - get information of track position and direction 
        else 
        {
            // First get reconstructed direction from tracker
            const Vector& trackDirection = vertex->getDirection();
            const Point&  trackPosition  = vertex->getPosition();
      
            double thetarec=-1;//reconstructed theta
            double phirec=0;  //reconstructed phi
            double px=0;      //projected position on (0,0,0 plane) 
            double py=0;
        
            thetarec=trackDirection[2];
            phirec=atan2(trackDirection[1],trackDirection[0]);
        
            //find projected positions on the 0,0,0 plane
            px = trackPosition[0]+sqrt((-trackPosition[2]/thetarec)*(-trackPosition[2]/thetarec)-trackPosition[2]*trackPosition[2])*cos(phirec) ; 
            py = trackPosition[1]+sqrt((-trackPosition[2]/thetarec)*(-trackPosition[2]/thetarec)-trackPosition[2]*trackPosition[2])*sin(phirec) ; 
            px = px+(m_k0-m_k1*thetarec+m_k2*thetarec*thetarec)*cos(phirec) ;
            py = py+(m_k0-m_k1*thetarec+m_k2*thetarec*thetarec)*sin(phirec) ;

                //log << MSG::INFO << "track direction = " << thetarec << endreq;

            //for now valid up to 26 degrees.. 
            if (-thetarec < 0.898 ) 
            {
                llStatus = -3. ;
                //return sc ;
            }
            else 
            {
                //make cuts for the 0-26 degree limit. Will later make the (theta,phi) dependence available
    
                if ( fabs(fabs(fabs(fabs(px)-375)-187.5)-187.5)>55 &&
                     fabs(fabs(fabs(fabs(py)-375)-187.5)-187.5)>55) 
                {
                    // Evaluation of energy using correlation method
                    // Coefficients fitted using GlastSim.
                    /*
                    double slope = trackDirection.z() ;
                    double p0 = -1.49 + 1.72 * slope ;
                    double p1 = 0.28 + 0.434 * slope ;
                    double p2 = -15.16 + 11.55 * slope ;
                    double p3 = 13.88 - 10.18 * slope ;
                    double lnE = log(eTotal/1000.);
                    double funcoef = (p0 + p1 * lnE )/(1+exp(-p2*(lnE - p3)));
                    */
          
                    double acoef = m_c0 - m_c1*thetarec + m_c2*thetarec*thetarec;
    
                    double lnE = log(eTotal/1000.) ;
                    double funcoef = (acoef + m_c3 * lnE ) ;
                    double E1 = eTotal + funcoef* (*cluster)[m_calNLayers-1].getEnergy();
                    double biascoef = m_b0 + m_b1*log(E1/1000);
                    double E2 = E1/biascoef;
                    funcoef = (acoef + m_c3 * log(E2/1000) ) ;
                    ///first iteration
                    double E3 = eTotal + funcoef* (*cluster)[m_calNLayers-1].getEnergy();
                    biascoef= m_b0 + m_b1*log(E3/1000);
                    double E4 = E3/biascoef;
                    funcoef = (acoef + m_c3 * log(E4/1000) );
                    ///second iteration.. probably overkill

                    params.setEnergy(E4);
    
                    // eCorrected = funcoef * cluster->getEneLayer()[getNLayers()-1] ;
                }
                else 
                {
                    llStatus = -4 ;
                }
            }
        }
    }

    // Whew! 
    // Ok, fill in the corrected information and exit
    
    corResult->setParams(params);
    corResult->insert(Event::CalCorEneValuePair("llStatus", llStatus));
    return corResult;
}
