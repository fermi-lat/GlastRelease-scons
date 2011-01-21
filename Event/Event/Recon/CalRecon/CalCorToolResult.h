#ifndef CalCorToolResult_H
#define CalCorToolResult_H

#include <vector>
#include <map>
#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "GaudiKernel/MsgStream.h"
#include "Event/Recon/CalRecon/CalParams.h"
#include "Event/RelTable/RelTable.h"

// static const CLID& CLID_CalCorToolResult = InterfaceID("CalCorToolResult", 1, 0);

/**
*  @class CalCorToolResult
*
*
*  @brief This defines the output from a given Energy Correction Tool 
*  
*  \author CalRecon Rewrite Group
*
* $Header$
*/

namespace Event 
{
/// Define a map to contain output of correction tools
/// Use a string for the key (to give flexibility), value to be contained a double
typedef std::map<std::string,double> CalCorEneValueMap;
typedef std::pair<std::string,double> CalCorEneValuePair;
    
class CalCorToolResult: public CalCorEneValueMap
{     
public:
        
    CalCorToolResult() : m_statusBits(0), m_chiSquare(0.)   {};
        
    virtual ~CalCorToolResult() {}

    /// Status word bits organized like:
    /// low:   |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
    ///         < volume info  >  <   hit type    > <    Energy Correction Status      >
    /// high:  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
    ///         <            Correction Algorithm Used             >  < more hit type >
    enum StatusBits {
    /// @param ZERO
                     ZERO         = 0x00000000, 
    /// @param VALIDPARAMS Bit set means particular energy correction is valid
                     VALIDPARAMS  = 0x00000001, 
// DC: redundant with correctionName
//    /// @param CALVALS    Bit set means Profile performed this correction
//                     CALVALS      = 0x00100000,  
//    /// @param PROFILE    Bit set means Profile performed this correction
//                     PROFILE      = 0x00200000,  
//    /// @param LASTLAYER  Bit set means Profile performed this correction
//                     LASTLAYER    = 0x00400000,  
//    /// @param RAWENERGY  Bit set means raw energy sum over all clusters
//                     RAWENERGY    = 0x20000000,  
//    /// @param GODOGSGO    Bit set means GODOGSGO performed this correction
//                     GODOGSGO     = 0x40000000
                                                  };  
// DC : better use checkStatusBit(VALIDPARAMS)
//    /// Answer quick questions based on status bits
//    inline bool validParams() const {return (m_statusBits & VALIDPARAMS)  == VALIDPARAMS;}

    /// Retrieve corrected information
    const std::string &  getCorrectionName() const {return m_correctionName;}
    const CalParams &    getParams()         const {return m_params;}
    double              getChiSquare()      const {return m_chiSquare;}

    /// 
    /// Start here the methods for setting the information 
    /// setCorrectionName for setting the corrector name
    /// @param name : the name (identifier) of the correction tool used
    inline void setCorrectionName(const std::string& name) {m_correctionName = name;}
    /// setParams for setting energy parameters
    /// @param params : the Cal energy parameters
    inline void setParams(const CalParams& params)         {m_params = params;}
    /// setChiSquare for setting the chisquare of the correction "fit"
    /// @param chiSquare the chiSquare resulting from this correction 
    inline void setChiSquare(double chiSquare)             {m_chiSquare = chiSquare;}

        /// Access individual status bits
    inline void setStatusBit( StatusBits bitToSet ) { m_statusBits |=  bitToSet ; }
    inline void clearStatusBit( StatusBits bitToClear ) { m_statusBits &= ~bitToClear ; }
    inline bool checkStatusBit( StatusBits bitToCheck ) const { return ((m_statusBits&bitToCheck)!=ZERO) ; }

    /// Access the status bits globally
    inline unsigned int getStatusBits() const { return m_statusBits ; }
    inline void setStatusBits( unsigned int statusBits ) { m_statusBits = statusBits ; }


private:
    std::string    m_correctionName;
    unsigned int   m_statusBits;
    double         m_chiSquare;
    CalParams      m_params;
};

//typedef for the Gaudi TDS Container
typedef std::vector<CalCorToolResult *> CalCorToolResultCol;
typedef CalCorToolResultCol::iterator CalCorToolResultColItr;
typedef CalCorToolResultCol::const_iterator CalCorToolResultColConItr;

// Define the relational table taking us back to CalCluster objects
class CalCluster;
typedef Event::RelTable<Event::CalCluster, Event::CalCorToolResult>               CalClusToCorToolHitTab;
typedef Event::Relation<Event::CalCluster, Event::CalCorToolResult>               CalClusToCorToolHitRel;
typedef ObjectList< Event::Relation<Event::CalCluster, Event::CalCorToolResult> > CalClustToCorToolHitTabList;

}

#endif        






