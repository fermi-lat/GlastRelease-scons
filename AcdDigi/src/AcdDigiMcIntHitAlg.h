#ifndef _AcdDigi_AcdDigiMcIntHitAlg_H
#define _AcdDigi_AcdDigiMcIntHitAlg_H 1

#include "GaudiKernel/Algorithm.h"

#include "AcdDigiUtil.h"
#include "idents/AcdId.h"

// to access an XML containing Digi parameters file
#include "xml/IFile.h"
#include "facilities/Util.h"


/** @class AcdDigiMcIntHitAlg
* @brief Algorithm to convert from hit data into digitization data 
* for the ACD.
* 
* This algorithm does the work of creating the AcdDigi objects and storing 
* them on the TDS, using the McIntegratingHit collection.
*
* @author Heather Kelly
* $Header$
*/

class AcdDigiMcIntHitAlg : public Algorithm {
    
public:
    
    AcdDigiMcIntHitAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
        
private:
    /// Read data from the input XML file
    void getParameters();
    
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
    
    /// standard deviation for gaussian noise for PHA, veto and CNO 
    /// discriminators
    double m_noise_std_dev_pha, m_noise_std_dev_veto, m_noise_std_dev_cno;
    
    /// full scale for PHA
    unsigned short m_full_scale;
    
    /// Global ratio of photoelectrons to mips
    unsigned short m_mean_pe_per_mip;
    
    /// number of MIPs that correspond to full scale PHA
    float m_mips_full_scale;
    
    /// MeV per MIP
    float m_mev_per_mip;
    
    /// jobOptions parameter denoting whether or not to perform auto 
    /// calibration to determine the number of MIPs for full scale PHA
    bool m_auto_calibrate;
    
    /// jobOptions parameter denoting whether or not to apply Poisson fluctuations
    bool m_apply_poisson;

    /// jobOptions parameter denoting whether or not to apply Gaussian noise 
    /// before determining PHA and discriminators
    bool m_apply_noise;
};

#endif