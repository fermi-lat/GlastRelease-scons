//                                                                       
// The CalDiagnosticData class contains the EM diagnostic CAL trigger primitives 
//                                                                       

#include "digiRootData/CalDiagnosticData.h"
#include <iostream>
// From Online EBF library
#include "CALdiagnostic.h"

ClassImp(CalDiagnosticData)

CalDiagnosticData::CalDiagnosticData() : m_datum(0)
{
}

CalDiagnosticData::CalDiagnosticData(UInt_t datum)
: m_datum(datum)
{
}

CalDiagnosticData::~CalDiagnosticData() {

}

void CalDiagnosticData::initialize(UInt_t datum)
{
    m_datum=datum;
}

void CalDiagnosticData::Clear(Option_t *option) {
    m_datum = 0;
}

void CalDiagnosticData::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << " Data Word: " << m_datum << endl;
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
