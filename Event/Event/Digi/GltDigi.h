#ifndef Event_GltDigi_H
#define Event_GltDigi_H 1

#include "GaudiKernel/ContainedObject.h"
#include <vector>

/*!
* \class GltDigi
* \author Johann Cohen-Tanugi
*
* \brief TDS class for Glt data
* TDS Trigger class for interface to filter algorithms
* 
* $Header$
*/

extern const CLID& CLID_GltDigi;

namespace Event {

class GltDigi : virtual public ContainedObject {
   
 public:
  GltDigi(){;}

  virtual ~GltDigi() {;}  

  inline std::vector<bool> getTkrThreeInRow(){return m_TKR_threeinRow;}
  inline std::vector<bool> getCAL_LO(){return m_CAL_LO;}
  inline std::vector<bool> getCAL_HI(){return m_CAL_HI;}

  inline void setTkrThreeInRow(std::vector<bool> value)
    {m_TKR_threeinRow = value;}

  inline void setCAL_LO(std::vector<bool> value)
    {m_CAL_LO = value;}

  inline void setCAL_HI(std::vector<bool> value)
    {m_CAL_HI = value;}

 private:
  std::vector<bool> m_TKR_threeinRow;
  std::vector<bool> m_CAL_LO;
  std::vector<bool> m_CAL_HI;
};
  
} //Namespace Event
#endif
