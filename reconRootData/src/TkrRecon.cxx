#include "reconRootData/TkrRecon.h"
#include "reconRootData/TkrTrack.h"
#include "reconRootData/TkrSiCluster.h"

ClassImp(TkrRecon)

//This is the default constructor (which does nothing but insure the pointers are zero)
TkrRecon::TkrRecon() 
{
    m_SiClusters = 0;
    m_TrackCands = 0;
    m_Tracks     = 0;
    m_Vertices   = 0;
}

//Destructor
TkrRecon::~TkrRecon() 
{
    Clean();

    // Delete the allocated memory
    if (m_SiClusters) {delete m_SiClusters;}
    if (m_TrackCands) {delete m_TrackCands;}
    if (m_Tracks)     {delete m_Tracks;}
    if (m_Vertices)   {delete m_Vertices;}

    // Zero the pointers (is this really necessary?)
    m_SiClusters = 0;
    m_TrackCands = 0;
    m_Tracks     = 0;
    m_Vertices   = 0;
}

// Clean out existing object arrays (e.g. between events)
void TkrRecon::Clean() 
{
    // First up clear the SiClusters
    if (m_SiClusters)
    {
        int nEntries = m_SiClusters->GetEntries();

        while(nEntries--) delete m_SiClusters->At(nEntries);

        m_SiClusters->Clear();
    }

    // Next are the pattern track candidates
    if (m_TrackCands)
    {
        int nEntries = m_TrackCands->GetEntries();

        while(nEntries--) delete m_TrackCands->At(nEntries);

        m_TrackCands->Clear();
    }

    // Now the fit tracks
    if (m_Tracks)
    {
        int nEntries = m_Tracks->GetEntries();

        while(nEntries--) delete m_Tracks->At(nEntries);

        m_Tracks->Clear();
    }

    // Now the vertices
    if (m_Vertices)
    {
        int nEntries = m_Vertices->GetEntries();

        while(nEntries--) delete m_Vertices->At(nEntries);

        m_Vertices->Clear();
    }
}

// Creates (i.e. allocates the memory for) the object arrays
void TkrRecon::Create() 
{
    if (!m_SiClusters) m_SiClusters = new TObjArray();
    if (!m_TrackCands) m_TrackCands = new TObjArray();
    if (!m_Tracks)     m_Tracks     = new TObjArray();
    if (!m_Vertices)   m_Vertices   = new TObjArray();
}
