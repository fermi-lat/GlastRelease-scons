#include "reconRootData/ReconHeader.h"

ClassImp(ReconHeader)

ReconHeader::ReconHeader() {
    m_flags = 0;
    m_eventId = 0;
    m_runId = 0;
}

int ReconHeader::Clear() {
    m_flags = 0;
    m_eventId = 0;
    m_runId = 0;
    return 0;
}

Bool_t ReconHeader::isValidTkrFlag(UInt_t val) {
    if ((val & (~tkrTestWord)) > 0) return kFALSE;
    return kTRUE;
}


Bool_t ReconHeader::setTkrFlag(UInt_t val) {

    if (!isValidTkrFlag(val)) return kFALSE;

    m_flags |= val;
    return kTRUE;
}


Bool_t ReconHeader::isValidCalFlag(UInt_t val) {
    if ((val & (~calTestWord)) > 0) return kFALSE;
    return kTRUE;
}


Bool_t ReconHeader::setCalFlag(UInt_t val) {

    if (!isValidCalFlag(val)) return kFALSE;

    m_flags |= val;
    return kTRUE;
}

