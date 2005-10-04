/** @file GltDigi.h
    @brief Definition and implementation of GltDigi

    $Header$


*/
#ifndef Event_GltDigi_H
#define Event_GltDigi_H 1

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

static const CLID& CLID_GltDigi = InterfaceID("GltDigi", 1, 0);

namespace Event {

class GltDigi : virtual public DataObject {
   
 public:
  GltDigi(){;}

  virtual ~GltDigi() {;}  

  inline std::vector<bool> getTkrThreeInRow()const{return m_TKR_threeinRow;}
  inline const std::vector<bool>& getCAL_LO()const{return m_CAL_LO;}
  inline const std::vector<bool>& getCAL_HI()const{return m_CAL_HI;}

  inline void setTkrThreeInRow(std::vector<bool> value)
    {m_TKR_threeinRow = value;}

  inline void setCAL_LO(std::vector<bool> value)
    {m_CAL_LO = value;}

  inline void setCAL_HI(std::vector<bool> value)
    {m_CAL_HI = value;}


    inline bool getCALLOtrigger()const { return triggerOR(m_CAL_LO);}
        

    inline bool getCALHItrigger()const {return triggerOR(m_CAL_HI);}

 private:
     inline bool triggerOR(const std::vector<bool>& bits)const{
        for( std::vector<bool>::const_iterator bit = bits.begin(); 
            bit!=bits.end(); ++ bit){
                if( *bit) { return true;}
            }
        return false;
     }

  std::vector<bool> m_TKR_threeinRow;
  std::vector<bool> m_CAL_LO;
  std::vector<bool> m_CAL_HI;
};
  
} //Namespace Event
#endif
