#ifndef CalEventEnergy_H
#define CalEventEnergy_H

#include <vector>
#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/MsgStream.h"

#include "Event/RelTable/RelTable.h"
#include "Event/Recon/CalRecon/CalCorToolResult.h"

//class Event::CalCluster;

static const CLID& CLID_CalEventEnergyCol = InterfaceID("CalEventEnergyCol", 1, 0);

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
    
class CalEventEnergy : public CalCorToolResultCol, virtual public ContainedObject
{     
public:
        
    CalEventEnergy() : m_statusBits(0), m_params() {}
        
    virtual ~CalEventEnergy() { clear() ; }

    void clear() ;

    /// Status word bits organized like:
    /// low:   |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
    ///         < volume info  >  <             hit type            > <  Recon Status  >
    /// high:  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
    ///         <            more bits to play with                 >  < more hit type >
    enum StatusBits {
    /// @param ZERO
                     ZERO = 0x00000000, 
    /// @param VALIDPARAMS Bit set means global params are valid
                     VALIDPARAMS  = 0x00000001, 
// DC: PASS_ONE and PASS_TWO can be avoided
//    /// @param PASS_ONE  Bit set means succesful first pass energy recon (no tracking)
//                     PASS_ONE  = 0x00000002, 
//    /// @param PASS_TWO  Bit set means successful second pass energy recon
//                     PASS_TWO  = 0x00000004,
//    /// @param GODOGSGO    Bit set means GODOGSGO is valid
//                     GODOGSGO     = 0x100000
    };  

//    /// Answer quick questions based on status bits
//    inline bool         validPassOne()  const {return (m_statusBits & PASS_ONE)  == PASS_ONE;}
//    inline bool         validPassTwo()  const {return (m_statusBits & PASS_TWO)  == PASS_TWO;}
//

    /// Find the most recent CalCorrToolResult from a given tool
    const CalCorToolResult * findLast( const std::string & correctionName ) const ;
    
    /// Access to "the" energy and parameters
    const CalParams&    getParams()     const {return m_params;}
    double              getEnergy()     const {return m_params.getEnergy();}
    const Point&        getCentroid()   const {return m_params.getCentroid();}
    const Vector&       getDirection()  const {return m_params.getAxis();}

    ///
    /// Set methods for this class
    /// @param energy the corrected energy
    inline void setParams(const CalParams& params)       {m_params = params;}

        /// Access individual status bits
    inline void setStatusBit( StatusBits bitToSet ) { m_statusBits |=  bitToSet ; }
    inline void clearStatusBit( StatusBits bitToClear ) { m_statusBits &= ~bitToClear ; }
    inline bool checkStatusBit( StatusBits bitToCheck ) const { return ((m_statusBits&bitToCheck)!=ZERO) ; }

    /// Access the status bits globally
    inline unsigned int getStatusBits() const { return m_statusBits ; }
    inline void setStatusBits( unsigned int statusBits ) { m_statusBits = statusBits ; }

        
private:
    unsigned int m_statusBits;
    CalParams    m_params;
};

//typedef for the Gaudi TDS Container
typedef ObjectVector<CalEventEnergy> CalEventEnergyCol ;
typedef CalEventEnergyCol::iterator CalEventEnergyColItr ;
typedef CalEventEnergyCol::const_iterator CalEventEnergyColConItr ;

}

#endif        






