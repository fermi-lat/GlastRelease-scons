// File and Version Information:
// $Header$
//
// Description: this utility singleton is used in various classes of G4Generator
// to register new McParticle objects, retrive the actual McParticle (i.e. the
// one that is actually creating hits in the sensitive detectors) and to finally
// save the McParticle hierarchy in the TDS
//      
// Author(s):
//      Iama Crazy

#include "DigiObjectManager.h"

// This is the singleton static pointer
DigiObjectManager* DigiObjectManager::m_pointer = 0;

// Define the default pool sizes for our objects
#define DIGIPOOLSIZE       100
#define DIGIMAXPOOLGROW    1000
#define DIAGPOOLSIZE       10
#define DIAGMAXPOOLGROW    1000
#define TEMPOOLSIZE        10
#define TEMMAXPOOLGROW     1000

DigiObjectManager::DigiObjectManager() : m_acdDigiPool(DIGIPOOLSIZE), 
                                         m_calDigiPool(DIGIPOOLSIZE),
                                         m_tkrDigiPool(DIGIPOOLSIZE),
                                         m_tkrDiagnosticDataPool(DIAGPOOLSIZE),
                                         m_calDiagnosticDataPool(DIAGPOOLSIZE),
                                         m_temPool(TEMPOOLSIZE)
{
    // Set pool iterators to the first item in each pool
    m_acdDigiPoolIdx           = m_acdDigiPool.begin();
    m_calDigiPoolIdx           = m_calDigiPool.begin();
    m_tkrDigiPoolIdx           = m_tkrDigiPool.begin();
    m_tkrDiagnosticDataPoolIdx = m_tkrDiagnosticDataPool.begin();
    m_calDiagnosticDataPoolIdx = m_calDiagnosticDataPool.begin();
    m_temPoolIdx               = m_temPool.begin();
}

DigiObjectManager* DigiObjectManager::getPointer()
{
  // Purpose and Method: standard singleton method to retrive the unique pointer
  if(m_pointer == 0) m_pointer = new DigiObjectManager();
  return m_pointer;
}

AcdDigi* DigiObjectManager::getNewAcdDigi()
{
    AcdDigi* acdDigi = 0;

    // If we have exceeded our pre-allocated list of AcdDigis then expand
    if (m_acdDigiPoolIdx == m_acdDigiPool.end())
    {
        // Add one more AcdDigi just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_acdDigiPoolIdx = m_acdDigiPool.insert(m_acdDigiPoolIdx, AcdDigi());

        // Expand the pool by some reasonable amount
        int newSize = std::max((int)(2 * m_acdDigiPool.size()), (int)DIGIMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_acdDigiPool.insert(m_acdDigiPool.end(), newSize, AcdDigi());
    }

    // Get the pointer to an available AcdDigi
    acdDigi = &*m_acdDigiPoolIdx++;

    return acdDigi;
}

CalDigi* DigiObjectManager::getNewCalDigi()
{
    CalDigi* calDigi = 0;

    // If we have exceeded our pre-allocated list of CalDigis then expand
    if (m_calDigiPoolIdx == m_calDigiPool.end())
    {
        // Add one more CalDigi just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_calDigiPoolIdx = m_calDigiPool.insert(m_calDigiPoolIdx, CalDigi());

        // Expand the pool by some reasonable amount
        int newSize = std::max((int)(2 * m_calDigiPool.size()), (int)DIGIMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_calDigiPool.insert(m_calDigiPool.end(), newSize, CalDigi());
    }

    // Get the pointer to an available CalDigi
    calDigi = &*m_calDigiPoolIdx++;

    return calDigi;
}

TkrDigi* DigiObjectManager::getNewTkrDigi()
{
    TkrDigi* tkrDigi = 0;

    // If we have exceeded our pre-allocated list of TkrDigis then expand
    if (m_tkrDigiPoolIdx == m_tkrDigiPool.end())
    {
        // Add one more TkrDigi just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_tkrDigiPoolIdx = m_tkrDigiPool.insert(m_tkrDigiPoolIdx, TkrDigi());

        // Expand the pool by some reasonable amount
        int newSize = std::max((int)(2 * m_tkrDigiPool.size()), (int)DIGIMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_tkrDigiPool.insert(m_tkrDigiPool.end(), newSize, TkrDigi());
    }

    // Get the pointer to an available TkrDigi
    tkrDigi = &*m_tkrDigiPoolIdx++;

    return tkrDigi;
}

TkrDiagnosticData* DigiObjectManager::getNewTkrDiagnosticData()
{
    TkrDiagnosticData* tkrDiagnosticData = 0;

    // If we have exceeded our pre-allocated list of TkrDiagnosticDatas then expand
    if (m_tkrDiagnosticDataPoolIdx == m_tkrDiagnosticDataPool.end())
    {
        // Add one more TkrDiagnosticData just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_tkrDiagnosticDataPoolIdx = m_tkrDiagnosticDataPool.insert(m_tkrDiagnosticDataPoolIdx, TkrDiagnosticData());

        // Expand the pool by some reasonable amount
        int newSize = std::max((int)(2 * m_tkrDiagnosticDataPool.size()), (int)DIAGMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_tkrDiagnosticDataPool.insert(m_tkrDiagnosticDataPool.end(), newSize, TkrDiagnosticData());
    }

    // Get the pointer to an available TkrDiagnosticData
    tkrDiagnosticData = &*m_tkrDiagnosticDataPoolIdx++;

    return tkrDiagnosticData;
}

CalDiagnosticData* DigiObjectManager::getNewCalDiagnosticData()
{
    CalDiagnosticData* calDiagnosticData = 0;

    // If we have exceeded our pre-allocated list of CalDiagnosticDatas then expand
    if (m_calDiagnosticDataPoolIdx == m_calDiagnosticDataPool.end())
    {
        // Add one more CalDiagnosticData just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_calDiagnosticDataPoolIdx = m_calDiagnosticDataPool.insert(m_calDiagnosticDataPoolIdx, CalDiagnosticData());

        // Expand the pool by some reasonable amount
        int newSize = std::max((int)(2 * m_calDiagnosticDataPool.size()), (int)DIAGMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_calDiagnosticDataPool.insert(m_calDiagnosticDataPool.end(), newSize, CalDiagnosticData());
    }

    // Get the pointer to an available CalDiagnosticData
    calDiagnosticData = &*m_calDiagnosticDataPoolIdx++;

    return calDiagnosticData;
}

Tem* DigiObjectManager::getNewTem()
{
    Tem* tem = 0;

    // If we have exceeded our pre-allocated list of Tems then expand
    if (m_temPoolIdx == m_temPool.end())
    {
        // Add one more Tem just before the end... this is really a 
        // trick to get a valid iterator at the end of the list
        m_temPoolIdx = m_temPool.insert(m_temPoolIdx, Tem());

        // Expand the pool by some reasonable amount
        int newSize = std::max((int)(2 * m_temPool.size()), (int)TEMMAXPOOLGROW);

        // For good measure expand the pool by the amount we just calculated
        m_temPool.insert(m_temPool.end(), newSize, Tem());
    }

    // Get the pointer to an available Tem
    tem = &*m_temPoolIdx++;

    return tem;
}

void DigiObjectManager::Delete()
{
    // Keep all of the pools down to a reasonable size...
    // Start with AcdDigi pool
    if (m_acdDigiPool.size() > 5*DIGIPOOLSIZE)
    {
        m_acdDigiPool.resize(5*DIGIPOOLSIZE);
    }

    // Now the CalDigi pool
    if (m_calDigiPool.size() > 5*DIGIPOOLSIZE)
    {
        m_calDigiPool.resize(5*DIGIPOOLSIZE);
    }

    // Now the TkrDigi pool
    if (m_tkrDigiPool.size() > 5*DIGIPOOLSIZE)
    {
        m_tkrDigiPool.resize(5*DIGIPOOLSIZE);
    }

    // Now the TkrDiagnosticData pool
    if (m_tkrDiagnosticDataPool.size() > 5*DIAGPOOLSIZE)
    {
        m_tkrDiagnosticDataPool.resize(5*DIAGPOOLSIZE);
    }

    // Now the CalDiagnosticData pool
    if (m_calDiagnosticDataPool.size() > 5*DIAGPOOLSIZE)
    {
        m_calDiagnosticDataPool.resize(5*DIAGPOOLSIZE);
    }

    // Now the Tem pool
    if (m_temPool.size() > 5*TEMPOOLSIZE)
    {
        m_temPool.resize(5*TEMPOOLSIZE);
    }

    // Ok, reset iterators to first element of our pools
    m_acdDigiPoolIdx       = m_acdDigiPool.begin();
    m_calDigiPoolIdx       = m_calDigiPool.begin();
    m_tkrDigiPoolIdx       = m_tkrDigiPool.begin();

    return;
}



