#include "reconRootData/CalRecon.h"
#include <iostream>

ClassImp(CalRecon)

CalRecon::CalRecon() {
   Clear();
}

CalRecon::~CalRecon() {
}

void CalRecon::Clear(Option_t *option) {
    m_xtalRecCol.clear();
    m_clusterCol.clear();
}

void CalRecon::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout << "# clusters: " << m_clusterCol.size()
        << "  # xtalRecData objects: " << m_xtalRecCol.size() << endl;
}