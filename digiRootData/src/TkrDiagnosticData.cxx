//                                                                       
// The TkrDiagnosticData class contains the information about the Level 1 Trigger count
// for an event.  
//                                                                       

#include "digiRootData/TkrDiagnostic.h"
#include <iostream>
// Online EBF Library
#include "TKRdiagnostic.h"

ClassImp(TkrDiagnosticData)

void TkrDiagnosticData::initialize(UInt_t datum)
{
    m_datum = datum;
}

void TkrDiagnosticData::Clear(Option_t *option) {
    m_datum = 0;
}

void TkrDiagnosticData::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Full data word " << m_datum << endl;
}

UInt_t TkrDiagnosticData::GTRC(Int_t gtrc) const {
    return TKRdiagnostic::GTRC(m_datum, gtrc);
}