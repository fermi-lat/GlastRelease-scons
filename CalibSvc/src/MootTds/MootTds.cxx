// $Header$

/** @file MootTds.cxx
    @author J. Bogart

    Non-inline implementation for moot-related objects in the tds.
*/

#include "CalibSvc/MootTds.h"

// copy constructors
MootBase::MootBase(const MootBase& other) : DataObject(other) {
  m_type = other.m_type;
  m_subtype = other.m_subtype;
}

MootParmCol::MootParmCol(MOOTSUBTYPE sub, unsigned key, unsigned sz)  :
   MootBase(MOOTTYPE_parm, sub), m_key(key) {

  if (sz) {
    m_ov.reserve(sz);

    // Now go ahead and create all the entries we'll need. They have
    // to be dummies for now
    for (unsigned ix = 0; ix < sz; ix++) {
      MootParm* dummy = new MootParm();
      m_ov.push_back(dummy);
    }
  }

}
MootParmCol::MootParmCol(const MootParmCol& other) : MootBase(other), 
                                                     m_key(other.m_key)
{
  m_ov = other.m_ov;
}

MootParmCol::~MootParmCol() {
  m_ov.clear();
}
