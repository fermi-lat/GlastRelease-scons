#include "reconRootData/TkrVertex.h"
#include <iostream>

ClassImp(TkrVertex)

TkrVertex::TkrVertex() 
{
    Clear();
}
    
TkrVertex::~TkrVertex() 
{
    Clear();
}

void TkrVertex::Clear(Option_t *option) 
{
    m_statusBits =  0;

    m_energy     =  0.;                     // energy associated with vertex
    m_position   = TVector3(0.,0.,0.);      // position of vertex
    m_direction  = TVector3(0.,0.,0.);      // direction of gamma causing pair conversion vertex

    m_chiSquare  =  0.;                     // Spacial chi-square for combining tracks
    m_quality    =  0.;                     // Vertex "Quality" derived from topology & chisq.
    m_arcLen1    =  0.;                     // Signed distance from head of track 1 to VTX
    m_arcLen2    =  0.;                     // Signed distance from head of track 1 to VTX
    m_doca       =  0.;                     // Distance between tracks at VTX location
    m_radlen     =  0.;                     // Integrated radiation lengths from end of track 1
    
    m_vtxID      = commonRootData::TkrId(); // Complete TkrId identifying the details of this vertex
		                                    // (This is the TkrId of the first hit after the vertex)
    m_params     = TkrTrackParams();        // Parameter structure for vertex (includes cov. matrix)

    m_tracks.Clear();
}

void TkrVertex::Print(Option_t *option) const {
    TObject::Print(option);
    using namespace std;
    cout << "Id: " << m_vtxID << endl;
    cout << "Pos: (" << m_position.X() << "," << m_position.Y()
        << "," << m_position.Z() << ")  Dir: (" 
        << m_direction.X() << "," << m_direction.Y() << ","
        << m_direction.Z() << ")" << endl;
}
