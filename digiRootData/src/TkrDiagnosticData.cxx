//                                                                       
// The TkrDiagnosticData class contains the information about the Level 1 Trigger count
// for an event.  
//                                                                       

#include "digiRootData/TkrDiagnosticData.h"
#include <iostream>
// Online EBF Library
#include "TKRdiagnostic.h"

ClassImp(TkrDiagnosticData)

void TkrDiagnosticData::initialize(UInt_t datum, UShort_t tower, UShort_t gtcc)
{
    m_datum = datum;
    m_tower = tower;  
    m_gtcc = gtcc;
}

void TkrDiagnosticData::Clear(Option_t *option) {
    m_datum = 0;
    m_tower = 0;
    m_gtcc = 0;
}

void TkrDiagnosticData::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Tower: " << m_tower << " GTCC: " << m_gtcc << " Full data word " << m_datum << endl;
}

UInt_t TkrDiagnosticData::GTRC(Int_t gtrc) const {
    return TKRdiagnostic::GTRC(m_datum, gtrc);
}
