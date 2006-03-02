//                                                                       
// The CalDiagnosticData class contains the EM diagnostic CAL trigger primitives 
//                                                                       

#include "digiRootData/CalDiagnosticData.h"
#include "Riostream.h"
// From Online EBF library
#include "CALdiagnostic.h"

ClassImp(CalDiagnosticData)

CalDiagnosticData::CalDiagnosticData() : m_datum(0),m_tower(0),m_layer(0)
{
}

CalDiagnosticData::CalDiagnosticData(UInt_t datum,UShort_t tower, UShort_t layer)
: m_datum(datum), m_tower(tower), m_layer(layer)
{
}

CalDiagnosticData::~CalDiagnosticData() {

}

void CalDiagnosticData::initialize(UInt_t datum, UShort_t tower, UShort_t layer)
{
    m_datum=datum;
    m_tower = tower;
    m_layer = layer;
}

    

void CalDiagnosticData::Clear(Option_t *option) {
    m_datum = 0;
    m_tower = 0;
    m_layer = 0;
}

void CalDiagnosticData::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Tower: " << m_tower << " layer: " << m_layer << " Data Word: " << m_datum << endl;
}

UInt_t CalDiagnosticData::high(UInt_t sign) const {
    return CALdiagnostic::high(m_datum, sign);
}

UInt_t CalDiagnosticData::low(UInt_t sign) const {
    return CALdiagnostic::low(m_datum, sign);
}

UInt_t CalDiagnosticData::logAccepts(UInt_t sign) const {
    return CALdiagnostic::logAccepts(m_datum, sign);
}
