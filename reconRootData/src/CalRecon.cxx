#include "reconRootData/CalRecon.h"
#include <iostream>

ClassImp(CalRecon)

CalRecon::CalRecon() {
   m_xtalRecCol = 0;
   m_clusterCol = 0;
}

CalRecon::~CalRecon() {
    if (m_clusterCol) {
        m_clusterCol->Delete();
        delete m_clusterCol;
        m_clusterCol = 0;
    }
    if (m_xtalRecCol) {
        m_xtalRecCol->Delete();
        delete m_xtalRecCol;
        m_xtalRecCol = 0;
    }
}

void CalRecon::initialize() {
    if (!m_clusterCol) m_clusterCol = new TObjArray();
    if (!m_xtalRecCol) m_xtalRecCol = new TObjArray();
}

void CalRecon::Clear(Option_t* /* option */) {
    if (m_xtalRecCol) m_xtalRecCol->Delete();
    if (m_clusterCol) m_clusterCol->Delete();
}

void CalRecon::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout << "# clusters: " << m_clusterCol->GetEntries()
        << "  # xtalRecData objects: " << m_xtalRecCol->GetEntries() << endl;
}
