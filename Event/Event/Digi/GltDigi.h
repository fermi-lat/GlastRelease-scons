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
#include <set>

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
    GltDigi(){;}

    virtual ~GltDigi() {;}  

    inline std::vector<bool> getTkrThreeInRow()const{return m_TKR_threeinRow;}

    inline void setTkrThreeInRow(std::vector<bool> value)
      {m_TKR_threeinRow = value;}

    /// return true if any Cal FLE trigger bit was set
    bool getCALLOtrigger() const {return (m_CAL_LO.size() > 0);}
    /// return true if any Cal FHE trigger bit was set
    bool getCALHItrigger() const {return (m_CAL_HI.size() > 0);}

    /** \brief return true if specified Cal FLE trigger bit was set
        \param xtalFaceId CalXtalId should include optional face information
     */
    bool getCALLOtrigger(const idents::CalXtalId xtalFaceId) const {
      return (m_CAL_LO.find(xtalFaceId) != m_CAL_LO.end());
    }

    /** \brief return true if specified Cal FHE trigger bit was set
        \param xtalFaceId CalXtalId should include optional face information
     */
    bool getCALHItrigger(const idents::CalXtalId xtalFaceId) const {
      return (m_CAL_HI.find(xtalFaceId) != m_CAL_HI.end());
    }


    /** \brief set specified Cal FLE trigger bit
        \param xtalFaceId CalXtalId should include optional face information
     */
    void setCALLOtrigger(const idents::CalXtalId xtalFaceId) {m_CAL_LO.insert(xtalFaceId);}

    /** \brief set specifiec Cal FHE trigger bit
        \param xtalFaceId CalXtalId should include optional face information
     */
    void setCALHItrigger(const idents::CalXtalId xtalFaceId) {m_CAL_HI.insert(xtalFaceId);}

    /// \brief return 16 bit trigger vector for FLE trigger, one bit per tower
    unsigned short getCALLETriggerVector() const;

    /// \brief return 16 bit trigger vector for FLE trigger, one bit per tower
    unsigned short getCALHETriggerVector() const;

  private:

    std::vector<bool> m_TKR_threeinRow;
    
    /** \brief contains set of all xtal faces with trigger = true
        \note CalXtalId should include optional face information
    */
    typedef std::set<idents::CalXtalId> CalTriggerSet;
    /// store all crystals with FLE trigger raised.
    CalTriggerSet m_CAL_LO;
    /// store all crystals with FHE trigger raised.
    CalTriggerSet m_CAL_HI;
  };
  
} //Namespace Event
#endif
