/** @file GltDigi.h
    @brief Definition and implementation of GltDigi

    $Header$


*/
#ifndef Event_GltDigi_H
#define Event_GltDigi_H 1

#include "idents/CalXtalId.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"
#include <vector>

/*!
 * \class GltDigi
 * \author Johann Cohen-Tanugi
 *
 * \brief TDS class for Glt data
 * TDS Trigger class for interface to filter algorithms
 * 
 */

static const CLID& CLID_GltDigi = InterfaceID("GltDigi", 2, 0);

namespace Event {

  class GltDigi : virtual public DataObject {
   
  public:
    GltDigi():
      m_CALLOTriggerVec(0),
      m_CALHITriggerVec(0)
    {;}

    virtual ~GltDigi() {;}  

    inline std::vector<bool> getTkrThreeInRow()const{return m_TKR_threeinRow;}

    inline void setTkrThreeInRow(std::vector<bool> value)
      {m_TKR_threeinRow = value;}

    /// return true if any Cal FLE trigger bit was set
    bool getCALLOtrigger() const {return (m_CALLOTriggerVec != 0);}
    /// return true if any Cal FHE trigger bit was set
    bool getCALHItrigger() const {return (m_CALHITriggerVec != 0);}


    /// store one bit for trigger in each of 16 cal modules (bits in
    /// tower number order)
    typedef unsigned short CalTriggerVec;
    
    /** \brief set specified Cal FLE trigger bit
        \param xtalFaceId CalXtalId should include optional face information
     */
    void setCALLOTriggerVec(const CalTriggerVec trigVec) {m_CALLOTriggerVec = trigVec;}

    /** \brief set specifiec Cal FHE trigger bit
        \param xtalFaceId CalXtalId should include optional face information
     */
    void setCALHITriggerVec(const CalTriggerVec trigVec) {m_CALHITriggerVec = trigVec;}

    /// \brief return 16 bit trigger vector for FLE trigger, one bit per tower
        CalTriggerVec getCALLOTriggerVec() const {return m_CALLOTriggerVec;}

    /// \brief return 16 bit trigger vector for FLE trigger, one bit per tower
        CalTriggerVec getCALHITriggerVec() const {return m_CALHITriggerVec;}

  private:

    std::vector<bool> m_TKR_threeinRow;

    /// store all crystals with FLE trigger raised.
    CalTriggerVec m_CALLOTriggerVec;
    /// store all crystals with FHE trigger raised.
    CalTriggerVec m_CALHITriggerVec;
  };
  
} //Namespace Event
#endif
