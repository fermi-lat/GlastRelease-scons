/** @file ObfStatus.h
* @author Tracy Usher
*
* $Header$

*/
#ifndef ObfStatus_H
#define ObfStatus_H 1

#include <map>
#include <set>
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/MsgStream.h"

#include "GaudiKernel/IInterface.h"

// Include these to make sure we have the status/veto bit definitions and masks
#include "EFC/GFC_V3_status.h"
#include "XFC/HFC_status.h"
#include "XFC/MFC_status.h"
#include "XFC/DFC_status.h"
#include "EDS/EDS_rsd.h"

// Include these to define the keys
#include "GFC_DB/GAMMA_DB_instance.h"
#include "XFC_DB/DGN_DB_instance.h"
#include "XFC_DB/HIP_DB_instance.h"
#include "XFC_DB/MIP_DB_instance.h"

static const CLID& CLID_ObfStatus = InterfaceID("ObfStatus", 1, 0);

namespace OnboardFilterTds {

/** 
* @class IObfStatus
*
* @brief Define an abtract class which will describe a base for the 
*        individual filter outputs. Currently, all filters output a 
*        32 bit status word but then also can include some number of
*        other items. This will allow us to store that info on the TDS
*
*/
class IObfStatus
{
public:
    // Virtual destructor
    virtual ~IObfStatus() {}

    // Retrieving the status word is common to all filters
    virtual unsigned int  getStatusWord()             const = 0;
    virtual int           getState()                  const = 0;
    virtual unsigned int  getVetoMask()               const = 0;
    virtual unsigned int  getVetoBit()                const = 0;

    // As is retrieving the id and status byte
    virtual unsigned char getFilterId()               const = 0;
    virtual unsigned char getFiltersb()               const = 0;

    // As is the word telling which prescaler affected the event
    virtual unsigned int  getPrescalerWord()          const = 0;

    // Output for to diagnostic stream 
    virtual std::ostream& fillStream(std::ostream& s) const = 0;
};


/** 
* @class ObfStatus
*
* @brief Describes the actual TDS object that contains the output status words 
*        from various filters run in the Onboard Filter environment
*
*/

class ObfFilterStatus : virtual public DataObject
{
public:

    // @brief Enumerate the keys for possible Filters we might have results for 
    // @brief The keys are meant to match the SchemaId's for each filter
    // @param GammaFilter - key for Gamma Filter results
    // @param CNOFilter   - key for CNO Filter results
    // @param MipFilter   - key for Mip Filter results
    enum FilterKeys {
        GammaFilter = GAMMA_DB_SCHEMA,
        HIPFilter   = HIP_DB_SCHEMA,
        MIPFilter   = MIP_DB_SCHEMA,
        DGNFilter   = DGN_DB_SCHEMA
    };

    // Standard constructor
    ObfFilterStatus();

    // Standard destructor
    virtual ~ObfFilterStatus();

    // Add results method
    void addFilterStatus(FilterKeys key, IObfStatus* status);

    // Return results method
    const IObfStatus* getFilterStatus(FilterKeys key) const;
        
    /// writes out the information of the cluster if msglevel is set to debug
    std::ostream& fillStream( std::ostream& s ) const;

private:
    typedef std::map<FilterKeys, IObfStatus*> ObfStatusMap;

    ObfStatusMap m_statusMap;
};

class ObfGammaStatus : virtual public IObfStatus
{
public:
    ObfGammaStatus(unsigned char id, unsigned int status, unsigned char sb, unsigned int prescaler, unsigned int energy) : 
      m_id(id), m_status(status), m_sb(sb), m_prescaler(prescaler), m_energy(energy) {}
    virtual ~ObfGammaStatus() {}

    // If msb is set below then event is to be vetoed
    unsigned int  getStatusWord()    const {return m_status;}
    int           getState()         const;

    unsigned int  getVetoMask()      const {return GFC_V3_STATUS_M_VETOES;}
    unsigned int  getVetoBit()       const {return GFC_V3_STATUS_M_VETOED;}

    unsigned char getFilterId()      const {return m_id;}
    unsigned char getFiltersb()      const {return m_sb;}
    unsigned int  getPrescalerWord() const {return m_prescaler;}

    unsigned int  getEnergyInLeus()  const;
    float         getEnergy()        const;
    unsigned int  getStageEnergy()   const {return m_energy;}
    unsigned int  getStage()         const {return (m_energy & GFC_V3_STAGE_M_STAGE) >> GFC_V3_STAGE_S_ENERGY;}
    
    std::ostream& fillStream(std::ostream& s) const
    {
        s << "Gamma Filter status:" << std::hex << m_status << std::dec << std::endl;
        return s;
    }
private:
    unsigned char m_id;
    unsigned int  m_status;
    unsigned char m_sb;
    unsigned int  m_prescaler;
    unsigned int  m_energy;
};

class ObfHipStatus : virtual public IObfStatus
{
public:
    ObfHipStatus(unsigned char id, unsigned int status, unsigned char sb, unsigned int prescaler) : 
      m_id(id), m_status(status), m_sb(sb), m_prescaler(prescaler) {}
    virtual ~ObfHipStatus() {}

     // If msb of below is set then event is to be vetoed
    unsigned int  getStatusWord()    const {return m_status;}
    int           getState()         const;

    unsigned int  getVetoMask()      const {return HFC_STATUS_M_VETO_DEF;}
    unsigned int  getVetoBit()       const {return HFC_STATUS_M_VETOED;}

    unsigned char getFilterId()      const {return m_id;}
    unsigned char getFiltersb()      const {return m_sb;}
    unsigned int  getPrescalerWord() const {return m_prescaler;}
    
    std::ostream& fillStream(std::ostream& s) const
    {
        s << "HIP Filter status:" << std::hex << m_status << std::dec << std::endl;
        return s;
    }
private:
    unsigned int  m_status;
    unsigned char m_id;
    unsigned char m_sb;
    unsigned int  m_prescaler;
};

class ObfMipStatus : virtual public IObfStatus
{
public:
    ObfMipStatus(unsigned char id, unsigned int status, unsigned char sb, unsigned int prescaler) : 
      m_id(id), m_status(status), m_sb(sb), m_prescaler(prescaler) {}
    virtual ~ObfMipStatus() {}

    unsigned int  getStatusWord()    const {return m_status;}
    int           getState()         const;

    unsigned int  getVetoMask()      const {return MFC_STATUS_M_VETO_DEF;}
    unsigned int  getVetoBit()       const {return MFC_STATUS_M_VETOED;}

    unsigned char getFilterId()      const {return m_id;}
    unsigned char getFiltersb()      const {return m_sb;}
    unsigned int  getPrescalerWord() const {return m_prescaler;}
    
    std::ostream& fillStream(std::ostream& s) const
    {
        s << "MIP Filter status:" << std::hex << m_status << std::dec << std::endl;
        return s;
    }
private:
    unsigned char m_id;
    unsigned int  m_status;
    unsigned char m_sb;
    unsigned int  m_prescaler;
};

class ObfDgnStatus : virtual public IObfStatus
{
public:
    ObfDgnStatus(unsigned char id, unsigned int status, unsigned char sb, unsigned int prescaler) : 
      m_id(id), m_status(status), m_sb(sb), m_prescaler(prescaler) {}
    virtual ~ObfDgnStatus() {}

    // If msb of below is set then event is to be vetoed
    unsigned int  getStatusWord()    const {return m_status;}
    int           getState()         const;

    unsigned int  getVetoMask()      const {return DFC_STATUS_M_VETO_DEF;}
    unsigned int  getVetoBit()       const {return DFC_STATUS_M_VETOED;}

    unsigned char getFilterId()      const {return m_id;}
    unsigned char getFiltersb()      const {return m_sb;}
    unsigned int  getPrescalerWord() const {return m_prescaler;}
    
    std::ostream& fillStream(std::ostream& s) const
    {
        s << "DGN Filter status:" << std::hex << m_status << std::dec << std::endl;
        return s;
    }
private:
    unsigned char m_id;
    unsigned int  m_status;
    unsigned char m_sb;
    unsigned int  m_prescaler;
};

} // end namespace

#endif // ObfStatus_H
