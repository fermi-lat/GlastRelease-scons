#include "reconRootData/TkrRecon.h"
#include "reconRootData/TkrTrack.h"
#include "reconRootData/TkrSiCluster.h"
#include <iostream>

ClassImp(TkrRecon)

TObjArray *TkrRecon::s_siClusterStaticCol=0;
TObjArray *TkrRecon::s_trackCandStaticCol=0;
TObjArray *TkrRecon::s_trackStaticCol=0;
TObjArray *TkrRecon::s_vertexStaticCol=0;


TkrRecon::TkrRecon() 
{
    if (!s_siClusterStaticCol) s_siClusterStaticCol = new TObjArray();
    m_siClusterCol = s_siClusterStaticCol;

    if (!s_trackCandStaticCol) s_trackCandStaticCol = new TObjArray();
    m_trackCandCol = s_trackCandStaticCol;

    if (!s_trackStaticCol) s_trackStaticCol = new TObjArray();
    m_trackCol = s_trackStaticCol;

    if (!s_vertexStaticCol) s_vertexStaticCol = new TObjArray();
    m_vertexCol = s_vertexStaticCol;

}

TkrRecon::~TkrRecon() 
{
    if (m_siClusterCol == s_siClusterStaticCol) s_siClusterStaticCol = 0;
    m_siClusterCol->Delete();
    delete m_siClusterCol;
    m_siClusterCol = 0;

    if (m_trackCandCol == s_trackCandStaticCol) s_trackCandStaticCol = 0;
    m_trackCandCol->Delete();
    delete m_trackCandCol;
    m_trackCandCol = 0;

    if (m_trackCol == s_trackStaticCol) s_trackStaticCol = 0;
    m_trackCol->Delete();
    delete m_trackCol;
    m_trackCol = 0;

    if (m_vertexCol == s_vertexStaticCol) s_vertexStaticCol = 0;
    m_vertexCol->Delete();
    delete m_vertexCol;
    m_vertexCol = 0;

}

void TkrRecon::Clear(Option_t *option) {
    m_siClusterCol->Delete();
    m_trackCandCol->Delete();
    m_trackCol->Delete();
    m_vertexCol->Delete();
}

void TkrRecon::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);

}