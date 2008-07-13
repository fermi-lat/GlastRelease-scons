#include "OnboardFilterTds/ObfFilterStatus.h"

// Use this to define the "state" of each filter
#include "enums/Lsf.h"

namespace OnboardFilterTds{
    
ObfFilterStatus::ObfFilterStatus()
{
    m_statusMap.clear();

    return;
}

// Standard destructor
ObfFilterStatus::~ObfFilterStatus()
{
    for (ObfStatusMap::iterator statIter = m_statusMap.begin(); statIter != m_statusMap.end(); statIter++)
    {
        IObfStatus* result = statIter->second;
        delete result;
    }
    m_statusMap.clear();

    return;
}

// Add results method
void ObfFilterStatus::addFilterStatus(FilterKeys key, IObfStatus* status)
{
    m_statusMap[key] = status;

    return;
}

// Return results method
const IObfStatus* ObfFilterStatus::getFilterStatus(FilterKeys key) const
{
    const IObfStatus* result = 0;

    ObfStatusMap::const_iterator statIter = m_statusMap.find(key);

    if (statIter != m_statusMap.end()) result = statIter->second;

    return result;
}

std::ostream& ObfFilterStatus::fillStream( std::ostream& s ) const
{
    // Purpose: writes out debug info
    // Inputs:  message stream
    // Outputs: data written to message stream
    for (ObfStatusMap::const_iterator statIter = m_statusMap.begin(); statIter != m_statusMap.end(); statIter++)
    {
        statIter->second->fillStream(s);
    }
/*
    s << " tray " << getTkrId().getTray() 
      << " face " << getTkrId().getBotTop()
      << " XY " << getTkrId().getView()
      << " pos (" << m_position.x() << ", " << m_position.y()
      << ", " << m_position.z() << ") "
      << " i0-if " << m_strip0 <<"-"<< m_stripf
      << " rawToT " << m_rawToT
      << " Mips " << m_ToT;
        //<< " status " << m_status; // std::hex << m_status << std::dec;
*/
    return s;
}

int  ObfGammaStatus::getState() const
{
    std::map<unsigned int, int> bitToEnumMap;

    bitToEnumMap[0x00] = enums::Lsf::PASSED;
    bitToEnumMap[0x20] = enums::Lsf::SUPPRESSED;
    bitToEnumMap[0x40] = enums::Lsf::VETOED;
    bitToEnumMap[0x60] = enums::Lsf::LEAKED;

    if (bitToEnumMap.find(m_sb) != bitToEnumMap.end()) return bitToEnumMap[m_sb];
    else return enums::Lsf::INVALID;

}

unsigned int ObfGammaStatus::getEnergyInLeus() const
{
    int energyBits = m_energy & GFC_STAGE_M_ENERGY;

    // Drag the sign bit up and down
    energyBits <<= 8;
    energyBits >>= 8;

    return energyBits;
}

float ObfGammaStatus::getEnergy() const
{
    int energyBits = m_energy & GFC_STAGE_M_ENERGY;

    // Drag the sign bit up and down
    energyBits <<= 8;
    energyBits >>= 8;

    float energy = static_cast<float>(energyBits) / 4.0;

    return energy;
}

int  ObfHipStatus::getState() const
{
    std::map<unsigned int, int> bitToEnumMap;

    bitToEnumMap[0x00] = enums::Lsf::PASSED;
    bitToEnumMap[0x20] = enums::Lsf::SUPPRESSED;
    bitToEnumMap[0x40] = enums::Lsf::VETOED;
    bitToEnumMap[0x60] = enums::Lsf::LEAKED;

    if (bitToEnumMap.find(m_sb) != bitToEnumMap.end()) return bitToEnumMap[m_sb];
    else return enums::Lsf::INVALID;

}


int  ObfMipStatus::getState() const
{
    std::map<unsigned int, int> bitToEnumMap;

    bitToEnumMap[0x00] = enums::Lsf::PASSED;
    bitToEnumMap[0x20] = enums::Lsf::SUPPRESSED;
    bitToEnumMap[0x40] = enums::Lsf::VETOED;
    bitToEnumMap[0x60] = enums::Lsf::LEAKED;

    if (bitToEnumMap.find(m_sb) != bitToEnumMap.end()) return bitToEnumMap[m_sb];
    else return enums::Lsf::INVALID;

}


int  ObfDgnStatus::getState() const
{
    std::map<unsigned int, int> bitToEnumMap;

    bitToEnumMap[0x00] = enums::Lsf::PASSED;
    bitToEnumMap[0x20] = enums::Lsf::SUPPRESSED;
    bitToEnumMap[0x40] = enums::Lsf::VETOED;
    bitToEnumMap[0x60] = enums::Lsf::LEAKED;

    if (bitToEnumMap.find(m_sb) != bitToEnumMap.end()) return bitToEnumMap[m_sb];
    else return enums::Lsf::INVALID;

}

} // namespace
