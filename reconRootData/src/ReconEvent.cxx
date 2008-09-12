#include "reconRootData/ReconEvent.h"
#include "Riostream.h"

ClassImp(ReconEvent)

ReconEvent::ReconEvent() {
    m_tkr = 0;
    m_cal = 0;
    m_acd = 0;
    m_adfRecon = 0;
    Clear();
}

ReconEvent::~ReconEvent() {

    if (m_cal) {
        delete m_cal;
        m_cal = 0;
    }

    if (m_tkr) {
        m_tkr->Clear();
        delete m_tkr;
        m_tkr = 0;
    }

    if (m_acd) {
        delete m_acd;
        m_acd = 0;
    }

    if (m_adfRecon) {
        delete m_adfRecon;
        m_adfRecon = 0;
    }

    Clear("ALL");

}

void ReconEvent::initialize(UInt_t eventId, UInt_t runId, TkrRecon *tkr, CalRecon *cal, AcdRecon *acd) {
    m_eventId = eventId;
    m_runId = runId;
    m_tkr = tkr;
    m_cal = cal;
    m_acd = acd;
}

void ReconEvent::initAdf(reconRootData::AdfRecon *adf) {
    m_adfRecon = adf;
}


void ReconEvent::Clear(Option_t* /* option */) {
    m_eventId = 0;
    m_runId = 0;
    if (m_cal) {
        delete m_cal;
        m_cal = 0;
    }
    if (m_tkr) {
        delete m_tkr;
        m_tkr = 0;
    }
    if (m_acd) {
        delete m_acd;
        m_acd = 0;
    }
    if (m_adfRecon) {
        delete m_adfRecon;
        m_adfRecon = 0;
    }
    m_eventFlags = 0;
    m_gleamEventFlags = 0;

    //ReconObjectManager::getPointer()->Delete(option);
}

void ReconEvent::Print(Option_t* /* option */) const {
    using namespace std;
    cout << "(Run, Event): (" << m_runId << ", " << m_eventId << ")" << endl;
    cout << "Fsw EventFlags: " << m_eventFlags << endl;
    cout << "Gleam EventFlags: " << m_gleamEventFlags << endl;
}


//======================================================
// For Unit Tests
//======================================================

void ReconEvent::Fake( Int_t /*ievent*/, Float_t /*randNum*/ ) {

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t ReconEvent::CompareInRange( const ReconEvent & ref, const std::string & name ) const {

    bool result = true ;

    //result = COMPARE_IN_RANGE(Energy) && result ;

    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;

}

