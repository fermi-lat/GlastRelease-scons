// File and Version Information:
//      $Header$
//
// Description:
//                                                
// Author:  Heather Kelly                       

#include "commonRootData/RelTable.h"
#include "Riostream.h"
#include "TCollection.h"  // Declares TIter

ClassImp(RelTable)

RelTable::RelTable() {
    m_table = 0;
    Clear();
}

RelTable::RelTable(const RelTable& copy) : TObject(copy) {
    Clear();
    TIter it(copy.m_table);
    Relation *rel = 0;
    while ((rel = (Relation*)it.Next())) {
        addRelation(rel);
    }
}

RelTable::~RelTable() {    
    Clear();
}


void RelTable::initialize(UInt_t nEvent, UInt_t nRun)
{
    Clear();
    m_table = new TObjArray();

    m_eventId = nEvent;
    m_runId   = nRun;
}

void RelTable::addRelation(Relation *rel) 
{
    if (!m_table) m_table = new TObjArray();
    m_table->Add(rel);
}

const Relation* RelTable::getRelByFirst(const TObject *key) const {
    TObjArrayIter it(m_table);
    Relation *rel = 0;
    while ((rel = (Relation*)it.Next())) {
        if (rel->getFirst() == key) return rel;
    }
    
    return 0;
}

void RelTable::Clear(Option_t *option) {
    TObject::Clear(option);
    if (m_table) 
    {
        m_table->SetOwner();
        m_table->Clear();
        delete m_table;
    }
    m_table = 0;
}

void RelTable::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    TIter it(m_table);
    Relation *rel = 0;
    while ((rel = (Relation*)it.Next())) {
        rel->Print();
    }
}


