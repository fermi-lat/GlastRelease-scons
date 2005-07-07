#include "reconRootData/CalRecon.h"
#include <iostream>

ClassImp(CalRecon)

CalRecon::CalRecon() 
{
    m_eventEnergyCol  = 0;
    m_xtalRecCol  = 0;
    m_clusterCol  = 0;
    m_mipTrackCol = 0;
}

CalRecon::~CalRecon() 
{
    if (m_eventEnergyCol) 
    {
        m_eventEnergyCol->Delete();
        delete m_eventEnergyCol;
        m_eventEnergyCol = 0;
    }
    
    if (m_clusterCol) 
    {
        m_clusterCol->Delete();
        delete m_clusterCol;
        m_clusterCol = 0;
    }
    
    if (m_xtalRecCol) 
    {
        m_xtalRecCol->Delete();
        delete m_xtalRecCol;
        m_xtalRecCol = 0;
    }
    
    if (m_mipTrackCol) 
    {
        m_mipTrackCol->Delete();
        delete m_mipTrackCol;
        m_mipTrackCol = 0;
    }
}

void CalRecon::initialize() 
{
    if (!m_eventEnergyCol)  m_eventEnergyCol  = new TObjArray();
    if (!m_clusterCol)  m_clusterCol  = new TObjArray();
    if (!m_xtalRecCol)  m_xtalRecCol  = new TObjArray();
    if (!m_mipTrackCol) m_mipTrackCol = new TObjArray();
}

void CalRecon::Clear(Option_t* /* option */) 
{
    if (m_eventEnergyCol)  m_eventEnergyCol->Delete();
    if (m_xtalRecCol)  m_xtalRecCol->Delete();
    if (m_clusterCol)  m_clusterCol->Delete();
    if (m_mipTrackCol) m_mipTrackCol->Delete();
}

void CalRecon::Print(Option_t *option) const 
{
    TObject::Print(option);
    std::cout
        << "# eventEnergies: " << m_eventEnergyCol->GetEntries()
        << " # mipTracks: " << m_mipTrackCol->GetEntries()
        << " # clusters: " << m_clusterCol->GetEntries()
        << " # xtalRecDatas: " << m_xtalRecCol->GetEntries() << std::endl;
}
