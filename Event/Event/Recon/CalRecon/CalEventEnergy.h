#ifndef CalEventEnergy_H
#define CalEventEnergy_H

#include <vector>
#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/MsgStream.h"

#include "Event/RelTable/RelTable.h"
#include "Event/Recon/CalRecon/CalCorToolResult.h"

class Event::CalCluster;

static const CLID& CLID_CalEventEnergy = InterfaceID("CalEventEnergy", 1, 0);

/**
*  @class CalEventEnergy
*
*
*  @brief This defines the top level CalRecon TDS object summarizing the 
*         Cal energy reconstruction. 
*         This is a new class and is intended to sit above CalCluster(s) in
*         the CalRecon hierarchy. 
*         The presumption is that there is but one of these objects, hence it
*         derives from Gaudi's DataObject. 
*  
*  \author CalRecon Rewrite Group
*
* $Header$
*/

namespace Event 
{
    
class CalEventEnergy: public CalCorToolResultCol //, public DataObject 
{     
public:
        
    CalEventEnergy() : m_statusBits(0), m_params() {}
        
    virtual ~CalEventEnergy() {}

    /// Status word bits organized like:
    /// low:   |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
    ///         < volume info  >  <             hit type            > <  Recon Status  >
    /// high:  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
    ///         <            more bits to play with                 >  < more hit type >
    /// @param PASS_ONE  Bit set means succesful first pass energy recon (no tracking)
    enum StatusBits {PASS_ONE  = 0x00000001, 
    /// @param PASS_TWO  Bit set means successful second pass energy recon
                     PASS_TWO  = 0x00000002,
    /// @param GODOGSGO    Bit set means GODOGSGO is valid
                     GODOGSGO     = 0x100000};  

    /// Access the status bits to determine details of the hit
    inline const unsigned int  getStatusBits() const {return m_statusBits;}

    /// Answer quick questions based on status bits
    inline const bool validPassOne() const {return (m_statusBits & PASS_ONE)  == PASS_ONE;}
    inline const bool validPassTwo() const {return (m_statusBits & PASS_TWO)  == PASS_TWO;}

    /// Access to "the" energy and parameters
    const CalParams  getParams()    const {return m_params;}
    const double     getEnergy()    const {return m_params.getEnergy();}
    const Point      getCentroid()  const {return m_params.getCentroid();}
    const Vector     getDirection() const {return m_params.getAxis();}

    ///
    /// Set methods for this class
    /// @param energy the corrected energy
    inline void setParams(const CalParams& params)    {m_params = params;}

    inline void setStatusBit(unsigned int bitToSet)   {m_statusBits |=  bitToSet;}
    inline void clearStatusBit(StatusBits bitToClear) {m_statusBits &= ~bitToClear;}
        
private:
    unsigned int m_statusBits;
    CalParams    m_params;
};

}

#endif	






