// $Header$
/** @class BadStrips
 *    Implementation of bad or hot strips TCDS representation
 */

#include "CalibData/Tkr/BadStrips.h"

namespace CalibData {
  BadStrips::Tower::Tower(bool allBad, int howBad, 
                          unsigned row, unsigned col) :
        m_allBad(allBad), m_howBad(howBad), m_row(row), m_col(col),
        m_uniplanes(0) {
    m_uniplanes = new std::vector<Uniplane>; 
  }

  BadStrips::Tower::Tower(const Tower& other) : m_allBad(other.m_allBad), 
    m_howBad(other.m_howBad), m_row(other.m_row), m_col(other.m_col),
    m_uniplanes(0)
  {
    m_uniplanes = new std::vector<Uniplane>;
    std::vector<Uniplane>::const_iterator iUni = (other.m_uniplanes)->begin();
    
    while (iUni != (other.m_uniplanes)->end() ) {
      m_uniplanes->push_back(*iUni);
      iUni++;
    }
  }

  BadStrips::Tower::~Tower() {

    //  Is this necessary, or does it happen already by virtue of the
    // clear?
    /*
      std::vector<Uniplane>::iterator iUni = m_uniplanes->begin();
    while (iUni != m_uniplanes->end() ) {
      delete iUni;
    }
    */

    m_uniplanes->clear();
    delete m_uniplanes;
  }

  BadStrips::Uniplane::Uniplane(bool allBad, int howBad, int tray, bool top,
                                const StripCol& strips) :
    m_allBad(allBad), m_howBad(howBad), m_tray(tray),
    m_top(top), m_badStrips(0) {
    m_badStrips = new StripCol(strips);
  }

  BadStrips::Uniplane::Uniplane(const Uniplane& other) : 
    m_allBad(other.m_allBad), m_howBad(other.m_howBad), m_tray(other.m_tray),
    m_top(other.m_top), m_badStrips(0)
  {
    m_badStrips = new StripCol(*(other.m_badStrips));
  }
    

  BadStrips::Uniplane::~Uniplane() {
    m_badStrips->clear();      // maybe not necessary but can't hurt
    delete m_badStrips;
  }

  BadStrips::BadStrips(eBadType bType, const ITime& since, const ITime& till, 
                       int serNo) :
    CalibBase(since, till, serNo), m_type(bType), m_towers(0)
  {
    m_towers = new std::vector<Tower>;
  }


  void BadStrips::update(CalibBase& other) {
    // The following dynamic_cast has got to work
    BadStrips& other1 = dynamic_cast<BadStrips& >(other);

    CalibBase::update(other1);
    m_type = other1.m_type;

    m_towers->clear();

    std::vector<Tower>::const_iterator iT = other1.m_towers->begin();

    while (iT != other1.m_towers->end() ) {
      m_towers->push_back(*iT);
      iT++;
    }
  }

  BadStrips::eBadType BadStrips::getBadType() const {return m_type;}

  eVisitorRet BadStrips::traverse(BadStripsVisitor *v) const {
    std::vector<Tower>::const_iterator iTower = m_towers->begin();
    eVisitorRet ret = DONE;
    while (iTower != m_towers->end() ) {
      if (iTower->m_allBad) {
        ret = v->badTower(iTower->m_row, iTower->m_col, iTower->m_howBad);
        if (ret != CONT) return ret;
      }
      //   If tower not all bad, loop over planes within towers
      else {
        std::vector<Uniplane>::const_iterator iUni = 
          (iTower->m_uniplanes)->begin();
        while (iUni != (iTower->m_uniplanes)->end() ) {
          ret = v->badPlane(iTower->m_row, iTower->m_col,
                            iUni->m_tray, iUni->m_top,
                            iUni->m_howBad, iUni->m_allBad,
                            *(iUni->m_badStrips));
          if (ret != CONT) return ret;
          iUni++;
        }
      }
      ++iTower;
    }
    // If got to here, traversed the entire data structure without
    // a murmur from client
    return DONE;
  }
  // No copy constructor for now.  If we ever have one, it should
  // be private.

  BadStrips::~BadStrips() {
    m_towers->clear();
    delete m_towers;
  }


  StatusCode BadStrips::addBadTower(bool allBad, int howBad, unsigned row,
                                unsigned col) {

    // Maybe first check to see if already have this tower??
    if (findTower(row, col)) {
      //  put out a message
      return StatusCode::FAILURE;
    }
    Tower newTower(allBad, howBad, row, col);
    m_towers->push_back(newTower);
    return StatusCode::SUCCESS;
  }

  StatusCode BadStrips::addBadPlane(unsigned short row, unsigned short col,
                                    unsigned int tray, bool top, int howBad,
                                    bool allBad, StripCol& badStrips) {
    // First check to see that
    //   we have specified tower in our list
    //   it's not marked as all bad
    Tower* pTower = findTower(row, col);
    if (!pTower) {
      // put out complaint
      return StatusCode::FAILURE;
    }
    if (pTower->m_allBad) {
      // put out complaint
      return StatusCode::FAILURE;
    }

    // If we already have this plane with same badness for this tower,
    // that also is an error (or we could merge the lists of strips..)
    std::vector<Uniplane>::const_iterator iUni=pTower->m_uniplanes->begin();

    while (iUni != pTower->m_uniplanes->end() ) {
      if ((iUni->m_tray == tray) && (iUni->m_top == top) && 
          (iUni->m_howBad == howBad)) {
        // complain
        return StatusCode::FAILURE;
      }
      iUni++;
    }

    // No conflict.  Go ahead and add the new uniplane
    Uniplane newUni(allBad, howBad, tray, top, badStrips);
    pTower->m_uniplanes->push_back(newUni);
    return StatusCode::SUCCESS;
  }

  BadStrips::Tower* BadStrips::findTower(unsigned row, unsigned col) {
    std::vector<Tower>::const_iterator iTower = m_towers->begin();

    while (iTower != m_towers->end() ) {
      if ((iTower->m_col == col) && (iTower->m_row == row)) {
        Tower* pTower = const_cast<Tower*>(iTower);
        return pTower;
      }
      iTower++;
    }
    return 0;
  }
  
}
