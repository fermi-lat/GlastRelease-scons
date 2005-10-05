/** @file GltDigi.h
    @brief Definition and implementation of GltDigi

    $Header$


*/
#ifndef Event_GltDigi_H
#define Event_GltDigi_H 1

#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"
#include <vector>
#include <map>

/*!
 * \class GltDigi
 * \author Johann Cohen-Tanugi
 *
 * \brief TDS class for Glt data
 * TDS Trigger class for interface to filter algorithms
 * 
 */

static const CLID& CLID_GltDigi = InterfaceID("GltDigi", 1, 0);

namespace Event {

  class GltDigi : virtual public DataObject {
   
  public:
    GltDigi(){;}

    virtual ~GltDigi() {;}  

    inline std::vector<bool> getTkrThreeInRow()const{return m_TKR_threeinRow;}

    /** \brief represents per xtal face Cal trigger info
        \note CalXtalId indeces should have optional face info set.
    */
    typedef std::map<idents::CalXtalId,bool> CalTriggerMap;
    inline const CalTriggerMap& getCAL_LO()const{return m_CAL_LO;}
    inline const CalTriggerMap& getCAL_HI()const{return m_CAL_HI;}

    inline void setTkrThreeInRow(std::vector<bool> value)
      {m_TKR_threeinRow = value;}

    inline void setCAL_LO(const CalTriggerMap &value)
      {m_CAL_LO = value;}

    inline void setCAL_HI(const CalTriggerMap &value)
      {m_CAL_HI = value;}


    inline bool getCALLOtrigger()const { return triggerOR(m_CAL_LO);}
        

    inline bool getCALHItrigger()const {return triggerOR(m_CAL_HI);}

  private:
    inline bool triggerOR(const CalTriggerMap& bits) const {
      for( CalTriggerMap::const_iterator bit = bits.begin(); 
           bit!=bits.end(); 
           ++ bit)
        {if( bit->second) { return true;}}

      return false;
    }

    std::vector<bool> m_TKR_threeinRow;
    CalTriggerMap m_CAL_LO;
    CalTriggerMap m_CAL_HI;
  };
  
} //Namespace Event
#endif
