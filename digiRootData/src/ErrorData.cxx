//                                                                       
// The ErrorData class contains the information about the error summary in
// on TEm 
//                                                                       

#include "digiRootData/ErrorData.h"
#include "Riostream.h"

ClassImp(ErrorData)


ErrorData::ErrorData() 
{
    Clear();
}

ErrorData::ErrorData(UShort_t c, UShort_t t, Bool_t phs, Bool_t tmo) {
    Clear();
    init(c, t, phs, tmo);
}

ErrorData::ErrorData(const ErrorData& t):TObject(t) {
    Clear();
    m_cal = t.m_cal;
    m_tkr = t.m_tkr;
    m_phs = t.m_phs;
    m_tmo = t.m_tmo;
}


void ErrorData::init(UShort_t c, UShort_t t, Bool_t phs, Bool_t tmo) {
    m_cal = c;
    m_tkr = t;
    m_phs = phs;
    m_tmo = tmo;
}


void ErrorData::Clear(Option_t *option) {
    m_cal = 0;
    m_tkr = 0;
    m_phs = kFALSE; 
    m_tmo = kFALSE;
}

void ErrorData::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Cal: " << m_cal << " Tkr: " << m_tkr << endl;
    cout << "Phased Error: " << m_phs << " TimeOut: " << m_tmo << endl;
   
}


