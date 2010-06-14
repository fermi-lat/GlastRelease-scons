#include "reconRootData/CalLogEne.h"

ClassImp(CalLogEne)

CalLogEne::CalLogEne() {
    m_id = 0;
}


CalLogEne::CalLogEne (UInt_t layer, UInt_t column, LogID::CALAxes xy) {
    m_id = new LogID();
    m_id->setColumn(column);
    m_id->setLayer(layer);
    m_id->setXY(xy);
}

CalLogEne::~CalLogEne() {
    Clean();
}


void CalLogEne::Clean() {
    if (m_id) {
        delete m_id;
        m_id = 0;
    }
}






