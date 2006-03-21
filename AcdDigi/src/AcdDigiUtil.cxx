#define AcdDigi_AcdDigiUtil_CPP 

// File and Version Information:
// $Header$
// Description
// Some utility methods helpful for performing the ACD digitization.

#include "AcdDigiUtil.h"

#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandGauss.h"

// for min and floor functions
#include <algorithm>
#include <cmath>

#include <map>
#include <utility>

// to access an XML containing Digi parameters file
#include "xmlBase/IFile.h"
#include "facilities/Util.h"

xmlBase::IFile *AcdDigiUtil::m_ifile;
double AcdDigiUtil::m_noise_std_dev_pha;
double AcdDigiUtil::m_noise_std_dev_veto;
double AcdDigiUtil::m_noise_std_dev_cno;
unsigned short AcdDigiUtil::m_full_scale;
unsigned short AcdDigiUtil::m_mean_pe_per_mip;
unsigned short AcdDigiUtil::m_mean_pe_per_mip_ribbon;
double AcdDigiUtil::m_mips_full_scale;
double AcdDigiUtil::m_mev_per_mip;
std::map< unsigned int, std::pair<float, float> > AcdDigiUtil::m_pePerMipMap;
std::map< unsigned int, std::pair<double, double> > AcdDigiUtil::m_gainMap;


AcdDigiUtil::AcdDigiUtil() {    
    
    m_ifile = 0;
}


AcdDigiUtil::~AcdDigiUtil() {
    if (m_ifile) delete m_ifile;
}


void AcdDigiUtil::getParameters(const std::string &xmlFile) {
    
    // Purpose and Method:  Read in the parameters from the XML file using IFile
    
    m_ifile = new xmlBase::IFile(xmlFile.c_str());
        
    m_mean_pe_per_mip = m_ifile->getInt("global_constants", "mean_pe_per_mip");
    m_mean_pe_per_mip_ribbon = m_ifile->getInt("global_constants", "mean_pe_per_mip_ribbon");
    
    m_noise_std_dev_pha = m_ifile->getDouble("global_constants", "noise_std_dev_pha");
    m_noise_std_dev_veto = m_ifile->getDouble("global_constants", "noise_std_dev_veto");
    m_noise_std_dev_cno = m_ifile->getDouble("global_constants", "noise_std_dev_cno");
    
    m_full_scale = m_ifile->getInt("global_constants", "full_scale");
    
    m_mips_full_scale = m_ifile->getDouble("global_constants", "mips_full_scale");
    
    m_mev_per_mip = m_ifile->getDouble("global_constants", "mev_per_mip");
    
    return;
}

std::ostream& AcdDigiUtil::dumpMeanPePerPmt(std::ostream& s) const {

    s  << "Dump All Mean PE per PMT values read in from XML" << std::endl;
    std::map< unsigned int, std::pair<float, float> >::const_iterator mapIt;
    for (mapIt = m_pePerMipMap.begin(); mapIt != m_pePerMipMap.end(); mapIt++) {
        std::pair<float, float> pePerMip = (*mapIt).second;
        s << "ID: " << (*mapIt).first << " PMTA: " << pePerMip.first 
          << " PMTB: " << pePerMip.second << std::endl;
    }
    return s;
}


double AcdDigiUtil::convertMevToMips(double energy_mev) {
    return energy_mev / m_mev_per_mip;
}

void AcdDigiUtil::convertMipsToPhotoElectrons(const idents::AcdId &id, 
                                             double pmtA_mips, unsigned int &pmtA_pe, 
                                             double pmtB_mips, unsigned int &pmtB_pe) {
    // Purpose and Method:  Convert from MIPs to PhotoElectrons for both PMTs of an AcdId.
    //   First check to see if an individual pePerMip has been assigned for this PMT
    //   Check both in a map that stores the values already read in from the input XML file
    //   and checks the XML file itself.  If no pePerMip is provided for these PMTs, we use
    //   the global m_mean_pe_per_mip to compute the photoelectrons for each PMT.
    // Input: pmtA_mips, pmtB_mips - the number of MIPs for each PMT and the AcdId
    // Output:  pmtA_pe and pmtB_pe - the number of photoelectrons for each PMT
    
    // Check the map
    if (m_pePerMipMap.find(id.id()) != m_pePerMipMap.end()) {
        std::pair<float, float> pePerMip = m_pePerMipMap[id.id()];
        pmtA_pe = (unsigned int) floor(pmtA_mips * pePerMip.first);
        pmtB_pe = (unsigned int) floor(pmtB_mips * pePerMip.second);
        return;
    } 
    // Check the XML file
    std::string idStr;
    facilities::Util::itoa(id.id(), idStr);
    std::string pmtIdStr = idStr;
    if (m_ifile->contains("meanPePerMip", pmtIdStr.c_str())) {
        std::vector<double> pePerMipVec = m_ifile->getDoubleVector("meanPePerMip", pmtIdStr.c_str());
        m_pePerMipMap[id.id()] = std::make_pair(pePerMipVec[0], pePerMipVec[1]);
        pmtA_pe = (unsigned int) floor(pmtA_mips * pePerMipVec[0]);
        pmtB_pe = (unsigned int) floor(pmtB_mips * pePerMipVec[1]);
        return;
    }
    // Now use the global mean_pe_per_mip
    if (id.tile()) {
        pmtA_pe = (unsigned int) floor(pmtA_mips * m_mean_pe_per_mip );
        pmtB_pe = (unsigned int) floor(pmtB_mips * m_mean_pe_per_mip );
    } else if (id.ribbon()) {
        pmtA_pe = (unsigned int) floor(pmtA_mips * m_mean_pe_per_mip_ribbon );
        pmtB_pe = (unsigned int) floor(pmtB_mips * m_mean_pe_per_mip_ribbon );
    }
    return;
}

void AcdDigiUtil::convertPhotoElectronsToMips(const idents::AcdId &id, 
                                             unsigned int pmtA_pe, double &pmtA_mips, 
                                             unsigned int pmtB_pe, double &pmtB_mips) {
    // Purpose and Method:  Convert from PhotoElectrons to MIPs for both PMTs of an AcdId.
    //   First check to see if an individual pePerMip has been assigned for this PMT
    //   Check both in a map that stores the values already read in from the input XML file
    //   and checks the XML file itself.  If no pePerMip is provided for these PMTs, we use
    //   the global m_mean_pe_per_mip to compute the MIPs for each PMT.
    // Input:  pmtA_pe and pmtB_pe - the number of photoelectrons for each PMT
    // Output: pmtA_mips, pmtB_mips - the number of MIPs for each PMT and the AcdId
    
    // Check the map
    if (m_pePerMipMap.find(id.id()) != m_pePerMipMap.end()) {
        std::pair<float, float> pePerMip = m_pePerMipMap[id.id()];
        pmtA_mips = pmtA_pe / pePerMip.first;
        pmtB_mips = pmtB_pe / pePerMip.second;
        return;
    } 
    // Check the XML file
    std::string idStr;
    facilities::Util::itoa(id.id(), idStr);
    std::string pmtIdStr = idStr;
    if (m_ifile->contains("meanPePerMip", pmtIdStr.c_str())) {
        std::vector<double> pePerMipVec = m_ifile->getDoubleVector("meanPePerMip", pmtIdStr.c_str());
        m_pePerMipMap[id] = std::make_pair(pePerMipVec[0], pePerMipVec[1]);
        pmtA_mips = pmtA_pe / pePerMipVec[0];
        pmtB_mips = pmtB_pe / pePerMipVec[1];
        return;
    }
    // Now use the global mean_pe_per_mip
    if (id.tile()) {
        pmtA_mips = ((double) pmtA_pe) / ( m_mean_pe_per_mip);
        pmtB_mips = ((double) pmtB_pe) / ( m_mean_pe_per_mip);
    } else if (id.ribbon()) {
        pmtA_mips = ((double) pmtA_pe) / ( m_mean_pe_per_mip_ribbon);
        pmtB_mips = ((double) pmtB_pe) / ( m_mean_pe_per_mip_ribbon);
    }

    return;
    
}

long AcdDigiUtil::shootPoisson(double pmtPhotoElectrons) {
    // Pupose and Method:  Returns a value from a Poisson distribution,
    //   using the input number of photoelectrons as the mean
    // Input:
    // Output: 
    
    return CLHEP::RandPoisson::shoot(pmtPhotoElectrons);
}

double AcdDigiUtil::shootGaussian(double std_dev) {
    // Purpose and Method:  Returns a value from a gaussian distribution, with mean
    //   zero and standard deviation determined by the one input parameter.
    // Input:  Standard Deviation
    // Output:  A value obtained from the the Gaussian distribution.
    
    return CLHEP::RandGauss::shoot(0.0, std_dev);
}

void AcdDigiUtil::calcMipsToFullScale(const idents::AcdId& id, 
                                     double pmtA_mips, unsigned int pmtA_pe, double &pmtA_mipsToFullScale, 
                                     double pmtB_mips, unsigned int pmtB_pe, double &pmtB_mipsToFullScale) {
    
    // Check the map
    if (m_pePerMipMap.find(id.id()) != m_pePerMipMap.end()) {
        std::pair<float, float> pePerMip = m_pePerMipMap[id.id()];

        pmtA_mipsToFullScale = m_mips_full_scale * (( pePerMip.first) / (pmtA_pe / pmtA_mips) );
        pmtB_mipsToFullScale = m_mips_full_scale * (( pePerMip.second) / (pmtB_pe / pmtB_mips) );

        return;
    } 
    // Check the XML file
    std::string idStr;
    facilities::Util::itoa(id.id(), idStr);
    std::string pmtIdStr = idStr;
    if (m_ifile->contains("meanPePerMip", pmtIdStr.c_str())) {
        std::vector<double> pePerMipVec = m_ifile->getDoubleVector("meanPePerMip", pmtIdStr.c_str());
        m_pePerMipMap[id] = std::make_pair(pePerMipVec[0], pePerMipVec[1]);


        pmtA_mipsToFullScale = m_mips_full_scale * (( pePerMipVec[0] ) / (pmtA_pe / pmtA_mips) );
        pmtB_mipsToFullScale = m_mips_full_scale * (( pePerMipVec[1] ) / (pmtB_pe / pmtB_mips) );

        return;
    }

    // Or use the global mean_pe_per_mip

    if (id.tile()) {
        pmtA_mipsToFullScale = m_mips_full_scale * (( m_mean_pe_per_mip) / (pmtA_pe / pmtA_mips) );
        pmtB_mipsToFullScale = m_mips_full_scale * (( m_mean_pe_per_mip) / (pmtB_pe / pmtB_mips) );
    } else if (id.ribbon()) {
        pmtA_mipsToFullScale = m_mips_full_scale * (( m_mean_pe_per_mip_ribbon) / (pmtA_pe / pmtA_mips) );
        pmtB_mipsToFullScale = m_mips_full_scale * (( m_mean_pe_per_mip_ribbon) / (pmtB_pe / pmtB_mips) );

    }
    
    return;
}

void AcdDigiUtil::applyGains(const idents::AcdId &id, 
                            double &pmtA_mipsToFullScale, double &pmtB_mipsToFullScale) {
    
    // Check the map
    if (m_gainMap.find(id.id()) != m_gainMap.end()) {
        std::pair<double, double> gain = m_gainMap[id.id()];
        pmtA_mipsToFullScale = m_mips_full_scale * gain.first;
        pmtB_mipsToFullScale = m_mips_full_scale * gain.second;
        return;
    } 
    
    // Check the XML file
    std::string idStr;
    facilities::Util::itoa(id.id(), idStr);
    std::string gainIdStr = "Gain_" + idStr;
    if (m_ifile->contains("Gains", gainIdStr.c_str())) {
        std::vector<double> gainVec = m_ifile->getDoubleVector("gains", gainIdStr.c_str());
        m_gainMap[id.id()] = std::make_pair( gainVec[0], gainVec[1] );
        pmtA_mipsToFullScale = m_mips_full_scale * gainVec[0];
        pmtB_mipsToFullScale = m_mips_full_scale * gainVec[1];
        return;
    }
    
    // No Gains provided for this id, use global mipsToFullScale
    pmtA_mipsToFullScale = m_mips_full_scale;
    pmtB_mipsToFullScale = m_mips_full_scale;
    
    return;
}



unsigned short AcdDigiUtil::convertMipsToPha(double mips, double mipsToFullScale) {
    
    return static_cast<unsigned short>(std::min((double)floor(mips / mipsToFullScale * m_full_scale), (double)m_full_scale));
    
}

bool AcdDigiUtil::compareVolIds(const idents::VolumeIdentifier& tileId, 
                                const idents::VolumeIdentifier& screwVolId) {

    ///Purpose and Method:  Compare the volumeIdentfiers of a tile and a 
    /// AcdScrewSq to see if the AcdScrewSq occurs within the given tile

    if ((tileId[0] != 1) && (screwVolId[0] != 1) ) return false;
    unsigned int i;
    // compare the entries one by one to see if they are equal
    for (i = 0; i<5; i++) {
        if(tileId[i] != screwVolId[i]) return false;
    }
    return true;
}
