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
#include "EFC/GFC_status.h"
#include "XFC/HFC_status.h"
#include "XFC/MFC_status.h"

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
    virtual unsigned int getStatusHi()                const = 0;
    virtual unsigned int getStatusLo()                const = 0;
    virtual unsigned int getStatus32()                const = 0;

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
    // @param GammaFilter - key for Gamma Filter results
    // @param CNOFilter   - key for CNO Filter results
    // @param MipFilter   - key for Mip Filter results
    enum FilterKeys {
        GammaFilter = 0,
        CNOFilter   = 1,
        MipFilter   = 2
    };

    // Standard constructor
    ObfFilterStatus();

    // Standard destructor
    virtual ~ObfFilterStatus();

    // Add results method
    void addFilterStatus(FilterKeys key, IObfStatus* status);

    // Return results method
    const IObfStatus* getFilterStatus(FilterKeys key);
        
    /// writes out the information of the cluster if msglevel is set to debug
    std::ostream& fillStream( std::ostream& s ) const;

private:
    typedef std::map<FilterKeys, IObfStatus*> ObfStatusMap;

    ObfStatusMap m_statusMap;
};

class ObfGammaStatus : virtual public IObfStatus
{
public:
    ObfGammaStatus(unsigned int status) : m_status(status) {}
    virtual ~ObfGammaStatus() {}

    // This for backward compatibility...
    unsigned int getStatusHi() const {return m_status >> GFC_STATUS_V_TKR_LT_2_ELO;}
    unsigned int getStatusLo() const {return m_status & 0xFFFF;}
    // If msb is set below then event is to be vetoed
    unsigned int getStatus32() const {return m_status;}

    unsigned int getVetoMask() const {return GFC_STATUS_M_VETOES;}
    unsigned int getVetoBit()  const {return GFC_STATUS_M_VETOED;}
    
    std::ostream& fillStream(std::ostream& s) const
    {
        s << "Gamma Filter status:" << std::hex << m_status << std::dec << std::endl;
        return s;
    }
private:
    unsigned int m_status;
};

class ObfCNOStatus : virtual public IObfStatus
{
public:
    ObfCNOStatus(unsigned int status) : m_status(status) {}
    virtual ~ObfCNOStatus() {}

    unsigned int getStatusHi() const {return m_status >> 16;}
    unsigned int getStatusLo() const {return m_status & 0xFFFF;}
    // If msb of below is set then event is to be vetoed
    unsigned int getStatus32() const {return m_status;}

    unsigned int getVetoMask() const {return HFC_STATUS_M_VETO_DEF;}
    unsigned int getVetoBit()  const {return HFC_STATUS_M_VETOED;}
    
    std::ostream& fillStream(std::ostream& s) const
    {
        s << "CNO Filter status:" << std::hex << m_status << std::dec << std::endl;
        return s;
    }
private:
    unsigned int m_status;
};

class ObfMipStatus : virtual public IObfStatus
{
public:
    ObfMipStatus(unsigned int status) : m_status(status) {}
    virtual ~ObfMipStatus() {}

    unsigned int getStatusHi() const {return m_status >> 16;}
    unsigned int getStatusLo() const {return m_status & 0xFFFF;}
    unsigned int getStatus32() const {return m_status;}

    unsigned int getVetoMask() const {return MFC_STATUS_M_VETO_DEF;}
    unsigned int getVetoBit()  const {return MFC_STATUS_M_VETOED;}
    
    std::ostream& fillStream(std::ostream& s) const
    {
        s << "MIP Filter status:" << std::hex << m_status << std::dec << std::endl;
        return s;
    }
private:
    unsigned int m_status;
};

} // end namespace

#endif // ObfStatus_H
