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

#include "AcdTileList.h"

#include <map>

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

    void clear();

    /// Read data from the input XML file
    void getParameters();

    /// add noise to all PMTs that exist (currently for tiles only)
    void addNoise();

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

    AcdTileList m_tiles;
    
    /// standard deviation for gaussian noise for PHA, 
    /// veto and CNO discriminators
    double m_noise_std_dev_pha, m_noise_std_dev_veto, m_noise_std_dev_cno;
    
    /// full scale for PHA
    unsigned short m_full_scale;
    
    /// Global ratio of photoelectrons to mips
    unsigned short m_mean_pe_per_mip;
    
    /// number of MIPs tha correspond to full scale PHA
    double m_mips_full_scale;
    
    /// MeV per MIP
    double m_mev_per_mip;

    /// Distance (mm) cutoff for applying edge effects
    double m_max_edge_dist;
    
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

    /// Slope of the linear function used to estimate the edge effect
    double m_edge_slope;
    /// y-intercept of the linear function used to estimate the edge effect
    double m_edge_intercept;

    std::map<const idents::AcdId, double> m_energyDepMap;

    std::map<const idents::AcdId, double> m_pmtA_toFullScaleMap;
    std::map<const idents::AcdId, double> m_pmtA_phaMipsMap;
    std::map<const idents::AcdId, double> m_pmtA_vetoMipsMap;
    std::map<const idents::AcdId, double> m_pmtA_cnoMipsMap;

    std::map<const idents::AcdId, double> m_pmtB_toFullScaleMap;
    std::map<const idents::AcdId, double> m_pmtB_phaMipsMap;
    std::map<const idents::AcdId, double> m_pmtB_vetoMipsMap;
    std::map<const idents::AcdId, double> m_pmtB_cnoMipsMap;


};

#endif