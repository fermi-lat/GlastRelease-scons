#include "reconRootData/TkrHit.h"
#include "reconRootData/TkrSiCluster.h"
#include "reconRootData/TkrLocator.h"

ClassImp(TkrHit)

TkrHit::TkrHit() {
    // default constructor
    m_residual = 0.0f;    // residual of the silicon cluster
    m_chi2 = 0.0f;        // chi**2

    m_cluster = 0;
    m_locator = 0;
}


TkrHit::~TkrHit() {
    Clean();
}

void TkrHit::Clean() {
    // not deleting the TkrSiCluster, since it's
    // probably pointing into a list of TkrSiClusters

    // Should I be deleting this?  Maybe not...it could be
    // a TkrLocator in the list of TkrLocators used for a
    // particular track????
    if (m_locator) {
        delete m_locator;
        m_locator = 0;
    }

    m_residual = 0.0f;    // residual of the silicon cluster
    m_chi2 = 0.0f;        // chi**2
}

void TkrHit::Create() {
}
