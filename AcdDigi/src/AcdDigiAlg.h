#ifndef _AcdDigi_AcdDigiAlg_H
#define _AcdDigi_AcdDigiAlg_H 1

#include "GaudiKernel/Algorithm.h"

#include "AcdDigiUtil.h"

#include "idents/AcdId.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/MonteCarlo/McPositionHit.h"

// to access an XML containing Digi parameters file
#include "xml/IFile.h"
#include "facilities/Util.h"


/** @class AcdDigiAlg
* @brief Algorithm to convert from hit data stored as McPositionHits into digitization data 
* for the ACD.
* 
* @author Heather Kelly
*
* $Header$
*/

class AcdDigiAlg : public Algorithm {
    
public:
    
    AcdDigiAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
        
private:
    /// Read data from the input XML file
    void getParameters();

    /// Adjusts the deposited energy recorded in an ACD volume 
    /// based on the location of the hit
    double edgeEffect(const Event::McPositionHit *hit);
    
    /// Low discrim threshold which activates PHA
    double m_low_threshold_mips;  
    /// High discrim threshold for CNO
    double m_high_threshold_mips;
    /// Veto discrim threshold for nominal veto signal
    double m_veto_threshold_mips; 
    
    /// input XML file containing parameters for Digitization
    std::string	m_xmlFile;

    /// Access the methods in the AcdDigiUtil class
    AcdDigiUtil util;

    /// access to the Glast Detector Service to read in geometry 
    /// constants from XML files
    IGlastDetSvc *m_glastDetSvc;
    
    /// standard deviation for gaussian noise for PHA, 
    /// veto and CNO discriminators
    double m_noise_std_dev_pha, m_noise_std_dev_veto, m_noise_std_dev_cno;
    
    /// full scale for PHA
    unsigned short m_full_scale;
    
    /// Global ratio of photoelectrons to mips
    unsigned short m_mean_pe_per_mip;
    
    /// number of MIPs tha correspond to full scale PHA
    float m_mips_full_scale;
    
    /// MeV per MIP
    float m_mev_per_mip;
    
    /// JobOptions parameter denoting whether or not to perform auto 
    /// calibration to determine the number of MIPs for full scale PHA
    bool m_auto_calibrate;
    
    /// JobOptions parameter denoting whether or not to apply Poisson fluctuations
    bool m_apply_poisson;

    /// JobOptions parameter denoting whether or not to apply Gaussian noise 
    /// before determining PHA and discriminators
    bool m_apply_noise;

    /// JobOptions parameter denoting whether or not to apply edge effects
    /// according to the position of MC hits.
    bool m_edge_effect;
};

#endif