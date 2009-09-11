//
// The Event class is the top level class for the ASCII to ROOT
// conversion. It contains all the information about a single event.
//

#include "digiRootData/EventSummaryData.h"
#include "Riostream.h"
#include "eventSummary.h"

using namespace std;

ClassImp(EventSummaryData)

EventSummaryData::EventSummaryData() {
    Clear();
}

EventSummaryData::~EventSummaryData() {
    Clear();
}

void EventSummaryData::initialize(UInt_t summary) {
    m_summary = summary;
}

void EventSummaryData::Clear(Option_t* option) {

    TObject::Clear(option);
    m_summary = 0;
    m_flags = 0;
    UInt_t i;
    for (i = 0; i<16; i++) { 
        m_temLen[i] = 0;
        m_diagLen[i] = 0;
        m_errLen[i] = 0;
    }
    m_otherContribLen[GEM] = 0;
    m_otherContribLen[AEM] = 0;
    m_otherContribLen[OSW] = 0;
    m_eventSizeInBytes = 0;
}

void EventSummaryData::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout<<"Event Size: " << m_eventSizeInBytes << endl;
    cout << "Tem Lengths: " << endl;
    UInt_t i;
    for (i = 0; i<16; i++)
        cout << "Tem: " << i << " Len: " << m_temLen[i] << endl;
        cout << "Diag: " << i << " Len: " << m_diagLen[i] << endl;
        cout << "Error: " << i << " Len: " << m_errLen[i] << endl;
    cout << "GEM Len: " << m_otherContribLen[GEM] << " AEM Len: " << 
         m_otherContribLen[AEM] << endl;
    cout << "OSW Len: " << m_otherContribLen[OSW] << " ERR Len: " << endl; 
}

UInt_t EventSummaryData::calStrobe() const {return EventSummary::calStrobe(m_summary); }

UInt_t EventSummaryData::tag() const { return EventSummary::tag(m_summary); }

UInt_t EventSummaryData::TACK() const { return EventSummary::TACK(m_summary); }

UInt_t EventSummaryData::readout4() const { return EventSummary::readout4(m_summary); }

UInt_t EventSummaryData::zeroSuppress() const { return EventSummary::zeroSuppress(m_summary); }

UInt_t EventSummaryData::marker() const { return EventSummary::marker(m_summary); }

UInt_t EventSummaryData::error() const { return EventSummary::error(m_summary); }

UInt_t EventSummaryData::diagnostic() const { return EventSummary::diagnostic(m_summary); }

UInt_t EventSummaryData::eventNumber() const { return EventSummary::eventNumber(m_summary); }

UInt_t EventSummaryData::getTrgParityError() const { return EventSummary::trgParityError(m_summary); }



