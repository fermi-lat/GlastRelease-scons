// File and Version Information:
//      $Header$
//
// Description:
//                                                
// Author:  Heather Kelly                       

#include "commonRootData/Relation.h"
#include <iostream>
#include "TCollection.h"  // Declares TIter

ClassImp(Relation)

Relation::Relation() {
    Clear();
}

Relation::Relation(TRef key, const TRefArray& valueCol) {
    m_key = key;
    UInt_t irel;
    for (irel = 0; irel < valueCol.GetEntries(); irel++) {
        m_valueCol.Add(valueCol.At(irel));
    }
}

Relation::~Relation() {    
    Clear();
}


void Relation::initialize(TRef key, const TRefArray& valueCol) {
    m_key = key;
    UInt_t irel;
    for (irel = 0; irel < valueCol.GetEntries(); irel++) {
        m_valueCol.Add(valueCol.At(irel));
}}

void Relation::Clear(Option_t *option) {
    m_key = 0;
    m_valueCol.Clear();
}

void Relation::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout << "Key: ";
    m_key.GetObject()->Print();
    cout << endl;
    
    UInt_t irel;
    for (irel = 0; irel < m_valueCol.GetEntries(); irel++) {
        (m_valueCol.At(irel))->Print();
        cout << endl;
    }
}


