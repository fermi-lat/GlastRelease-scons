#include "reconRootData/TkrRecon.h"

ClassImp(TkrRecon)

TkrRecon::TkrRecon() {

}

TkrRecon::~TkrRecon() 
{
    Clear();


}

void TkrRecon::Clear(Option_t *option) {
    m_siClusterCol.clear();
    m_trackCandCol.clear();
   // m_trackCol.clear();
 //   m_vertexCol.clear();
}

void TkrRecon::Print(Option_t *option) const {
}