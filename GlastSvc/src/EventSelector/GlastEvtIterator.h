// $Header$
//====================================================================
//  GlastEvtIterator.h
//--------------------------------------------------------------------
//
//
//
//====================================================================

#ifndef GLASTCNV_GLASTEVENTSELECTORITERATOR_H
#define GLASTCNV_GLASTEVENTSELECTORITERATOR_H 1

// Include files
#include "Gaudi/Interfaces/IEvtSelector.h"
#include "Gaudi/Interfaces/ISvcLocator.h"
#include "Gaudi/Kernel/Service.h"

#include <list>
#include <string>

// Forward declarations

class GlastEvtIterator: public IEvtSelector::Iterator {
public:
  friend class EventSelector;
  GlastEvtIterator() {}
  GlastEvtIterator( const IEvtSelector* glastEvtSel): m_glastEvtSel(glastEvtSel),
    m_recId(0) { }
  GlastEvtIterator( const IEvtSelector* glastEvtSel, int recId, std::list<std::string>::const_iterator inpIt)
    : m_glastEvtSel(glastEvtSel){
	  m_inputDataIt =  inpIt;
	  m_recId = recId;
  }
  virtual ~GlastEvtIterator() { }
  virtual IOpaqueAddress* operator*() const  { return m_glastEvtSel->reference(*this);}
  virtual IOpaqueAddress* operator->() const { return m_glastEvtSel->reference(*this);}
  virtual IEvtSelector::Iterator& operator++()    { return m_glastEvtSel->next(*this);}
  virtual IEvtSelector::Iterator& operator--()    { return m_glastEvtSel->previous(*this);}
  virtual IEvtSelector::Iterator& operator++(int) { return m_glastEvtSel->next(*this);}
  virtual IEvtSelector::Iterator& operator--(int) { return m_glastEvtSel->previous(*this);}
  virtual bool operator==(const IEvtSelector::Iterator& it) const {
    const GlastEvtIterator* glastIt = dynamic_cast<const GlastEvtIterator*>(&it);
    return((m_recId == glastIt->m_recId));
  }
  virtual bool operator!=(const IEvtSelector::Iterator& it) const {
    const GlastEvtIterator* glastIt = dynamic_cast<const GlastEvtIterator*>(&it);
    return((m_recId != glastIt->m_recId));
  }

private:
  const IEvtSelector*   m_glastEvtSel;
  int                   m_recId;
  std::list<std::string>::const_iterator  m_inputDataIt;

};

#endif  // GLASTCNV_GLASTEVENTSELECTORITERATOR_H



















