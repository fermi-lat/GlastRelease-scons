#ifndef _AcdDigi_AcdDigiUtil_H
#define _AcdDigi_AcdDigiUtil_H 1

#include <map>
#include <utility>
#include "Event/Digi/AcdDigi.h"

#include "idents/AcdId.h"

// to access an XML containing Digi parameters file
#include "facilities/Util.h"

//namespace xml { class IFile; }
#include "xml/IFile.h"

/** @class AcdDigiUtil
* @brief Utility class that defines the methods used for ACD digitization
* for the ACD.
* 
*
* $Header$
*/

class AcdDigiUtil  {
    
public:
    
    AcdDigiUtil(); 

    ~AcdDigiUtil();

    /// Read data from the input XML file
    void getParameters(const std::string &xmlFile);

    /// Convert from energy (MeV) to MIPs
    static double convertMevToMips(double energy_mev);

    /// Converts MIPs to PhotoElectrons
    /// @param id an AcdId identifer
    /// used to retrieve specific conversion parameters if available
    /// @param pmtA_mips input number of MIPs detected by PMT A
    /// @param pmtA_pe output number of PEs detected by PMT A
    /// @param pmtB_mips input number of MIPs detected by PMT B
    /// @param pmtB_pe output number of PEs detected by PMT B
    static void convertMipsToPhotoElectrons(const idents::AcdId &id, 
        float pmtA_mips, unsigned int &pmtA_pe,
        float pmtB_mips, unsigned int &pmtB_pe);

    /// Converts PhotoElectrons to MIPs
    /// @param id an AcdId identifer
    /// used to retrieve specific conversion parameters if available
    /// @param pmtA_pe input number of PEs detected by PMT A
    /// @param pmtA_mips output number of MIPs detected by PMT A
    /// @param pmtB_pe input number of PEs detected by PMT B
    /// @param pmtB_mips output number of MIPs detected by PMT B
    static void convertPhotoElectronsToMips(const idents::AcdId &id, 
        unsigned int pmtA_pe, float &pmtA_mips,
        unsigned int pmtB_pe, float &pmtB_mips);

    /// Determine MIPS to Full Scale conversion for both PMTs for a particular AcdId
    static void calcMipsToFullScale(const idents::AcdId&, 
        float pmtA_mips, unsigned int pmtA_pe, float &mipsToFullScaleA, 
        float pmtB_mips, unsigned int pmtB_pe, float& mipsToFullScaleB );

    /// calculate mipsToFullScale based on the specific parameters provided for AcdIds
    /// @param id an AcdId identifer
    /// @param mipsToFullScaleA the calculated parameter relating MIPs to Full Scale for PMT A
    /// @param mipsToFullScaleB the calculated parameter relating MIPs to Full Scale for PMT B
    static void applyGains(const idents::AcdId &id, 
        float &mipsToFullScaleA, float &mipsToFullScaleB);
    
    /// Converts MIPs to a PHA value - if the number of MIPs is offscale, returns
    /// the fullscale value
    static unsigned short convertMipsToPha(float mips, float mipsToFullScale);

    /// Returns a value sampled from a Poisson distribution
    /// @param pmtPhotoElectrons is the mean of the Poisson distribution
    static long calcPoisson(long pmtPhotoElectrons);

    /// Returns a value sampled from a Gaussian distribution
    /// @param std_dev Standard Deviation to be  used when sampling
    static float calcGaussianNoise(float std_dev);
    
private:
            
    static xml::IFile *m_ifile;
    
    /// standard deviation for gaussian noise for PHA, veto and CNO discriminators
    static double m_noise_std_dev_pha, m_noise_std_dev_veto, m_noise_std_dev_cno;
    
    /// full scale for PHA
    static unsigned short m_full_scale;
    
    /// Global ratio of photoelectrons to mips
    static unsigned short m_mean_pe_per_mip;
    
    /// number of MIPs tha correspond to full scale PHA
    static float m_mips_full_scale;
    
    /// MeV per MIP
    static float m_mev_per_mip;
     
    /// store AcdId specific number of photoelectrons per mip as they are read in
    static std::map< unsigned int, std::pair<int, int> > m_pePerMipMap;
    /// store AcdId specific gain values as they are read in
    static std::map< unsigned int, std::pair<float, float> > m_gainMap;
};

#endif