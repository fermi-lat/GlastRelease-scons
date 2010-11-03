// File and Version Information:
//      $Header$
//
// Description:
//                                                
// Author:  Heather Kelly                       

#include "commonRootData/Relation.h"
#include "Riostream.h"
#include "TCollection.h"  // Declares TIter

ClassImp(Relation)

Relation::Relation() {
    Clear();
}

Relation::Relation(TRef& first, TRef& second, const TObjArray& infos) 
{
    initialize(first, second, infos);
}

Relation::~Relation() {    
    Clear();
}


void Relation::initialize(TRef& first, TRef& second, const TObjArray& infos) 
{
    m_first  = first;
    m_second = second;

    // We need to make our own copy that we "own"
    m_infos = *(TObjArray*)infos.Clone();
}

void Relation::Clear(Option_t *option) 
{
    TObject::Clear(option);
    m_first  = 0;
    m_second = 0;
    m_infos.SetOwner();
    m_infos.Clear();
}

void Relation::Print(Option_t *option) const 
{
    using namespace std;
    TObject::Print(option);
    cout << "First: ";
    m_first.GetObject()->Print();
    cout << ", Second: ";
    m_second.GetObject()->Print();
    cout << endl;
    
    Int_t irel;
    for (irel = 0; irel < m_infos.GetEntries(); irel++) 
    {
        (m_infos.At(irel))->Print();
        cout << endl;
    }
}


