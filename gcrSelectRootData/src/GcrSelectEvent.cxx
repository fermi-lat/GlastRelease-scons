#include "gcrSelectRootData/GcrSelectEvent.h"
#include "Riostream.h"

ClassImp(GcrSelectEvent)

GcrSelectEvent::GcrSelectEvent() {
    m_gcrSelect = 0;
    Clear();
}

GcrSelectEvent::~GcrSelectEvent() {

    if (m_gcrSelect) {
        delete m_gcrSelect;
        m_gcrSelect = 0;
    }
    Clear();


}

void GcrSelectEvent::initialize(UInt_t eventId, UInt_t runId, GcrSelect* gcrSelect) {
    m_eventId = eventId;
    m_runId = runId;
    m_gcrSelect = gcrSelect;
}

void GcrSelectEvent::Clear(Option_t* /* option */) {
    m_eventId = 0;
    m_runId = 0;
    
    if (m_gcrSelect) {
        delete m_gcrSelect;
        m_gcrSelect = 0;
    }
    m_eventFlags = 0;
}

void GcrSelectEvent::Print(Option_t* /* option */) const {
    using namespace std;
    cout << "(Run, Event): (" << m_runId << ", " << m_eventId << ")" << endl;
    cout << "EventFlags: " << m_eventFlags << endl;
}


//======================================================
// For Unit Tests
//======================================================

/**void ReconEvent::Fake( Int_t , Float_t  ) {

}*/

// #define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

/**Bool_t ReconEvent::CompareInRange( const ReconEvent & ref, const std::string & name ) const {

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

}*/

