#include "reconRootData/CalRecon.h"
#include <iostream>

ClassImp(CalRecon)

CalRecon::CalRecon() {
}


CalRecon::~CalRecon() {
    m_clusterCol->Delete();
//    m_xtalRecCol->Delete();
    m_xtalRecCol.clear();
}

void CalRecon::Clear(Option_t *option) {

}

void CalRecon::Print(Option_t *option) {
    using namespace std;
    TObject::Print(option);
    cout << "# clusters: " << m_clusterCol->GetEntries() 
        << "  # xtalRecData objects: " << m_xtalRecCol.size() << endl;
}