#ifndef EBFEVTSELECTORITERATOR_H
#define EBFEVTSELECTORITERATOR_H 1

#include "GaudiKernel/IEvtSelector.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/Service.h"

#include <list>
#include <string>

/* @class LdfEventIterator
* @brief
*
* $Header$
*/
class LdfEvtIterator: public IEvtSelector::Iterator {
public:
  friend class LdfEventSelector;
  LdfEvtIterator() {}
  LdfEvtIterator( const IEvtSelector* glastEvtSel): m_glastEvtSel(glastEvtSel),
    m_recId(0) { }
  LdfEvtIterator(const IEvtSelector* glastEvtSel, int recId, 
                 std::list<std::string>::const_iterator inpIt)
    : m_glastEvtSel(glastEvtSel){
	  m_inputDataIt =  inpIt;
	  m_recId = recId;
          m_evtCount = 0;
  }
  virtual ~LdfEvtIterator() { }
  virtual IOpaqueAddress* operator*() const  
  { return m_glastEvtSel->reference(*this);}
  virtual IOpaqueAddress* operator->() const 
  { return m_glastEvtSel->reference(*this);}
  virtual IEvtSelector::Iterator& operator++()    
  { return m_glastEvtSel->next(*this);}
  virtual IEvtSelector::Iterator& operator--()    
  { return m_glastEvtSel->previous(*this);}
  virtual IEvtSelector::Iterator& operator++(int) 
  { return m_glastEvtSel->next(*this);}
  virtual IEvtSelector::Iterator& operator--(int) 
  { return m_glastEvtSel->previous(*this);}
  virtual bool operator==(const IEvtSelector::Iterator& it) const {
    const LdfEvtIterator* glastIt = dynamic_cast<const LdfEvtIterator*>(&it);
    return((m_recId == glastIt->m_recId));
  }
  virtual bool operator!=(const IEvtSelector::Iterator& it) const {
    const LdfEvtIterator* glastIt = dynamic_cast<const LdfEvtIterator*>(&it);
    return((m_recId != glastIt->m_recId));
  }
  void setFileNameIterator(std::list<std::string>::const_iterator inpIt) 
  { m_inputDataIt = inpIt; }
private:
  const IEvtSelector*   m_glastEvtSel;
  int                   m_recId;
  std::list<std::string>::const_iterator  m_inputDataIt;
  int                   m_evtCount;

};

#endif  // LdfEventSELECTORITERATOR_H



















