#include "reconRootData/TkrVertex.h"
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(TkrVertex)

TkrVertex::TkrVertex() 
{
    Clear();
}
    
TkrVertex::~TkrVertex() 
{
    Clear();
}

void TkrVertex::Clear(Option_t* /* option */) 
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


//======================================================
// For Unit Tests
//======================================================

void TkrVertex::Fake( Int_t ievent, UInt_t rank, Float_t randNum ) {

        Float_t f = Float_t(ievent) ;

        commonRootData::TkrId tkrId ;
        tkrId.Fake(ievent,rank,randNum) ;
        setTkrID(tkrId);
        
        Double_t qual = f;
        setQuality(qual);
        Double_t energy = f*randNum;
        setEnergy(energy);

        TkrTrackParams vtxPar ;
        vtxPar.Fake(ievent,rank,randNum) ;
        setParams(vtxPar);

        TVector3 pos(f, f*randNum, f*randNum*2.);
        setPosition(pos);
        TVector3 dir(randNum*2., randNum*3., randNum*4.);
        setDirection(dir);
        
        //addTrackId(rank+1);
        //addTrackId(rank+2);

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t TkrVertex::CompareInRange( const TkrVertex & ref, const std::string & name ) const {

    bool result = true ;

    result = COMPARE_IN_RANGE(Energy) && result ;
    result = COMPARE_IN_RANGE(Quality) && result ;
    result = COMPARE_IN_RANGE(TkrId) && result ;
    result = COMPARE_IN_RANGE(Position) && result ;
    result = COMPARE_IN_RANGE(Direction) && result ;
    result = COMPARE_IN_RANGE(VertexParams) && result ;

    //getNumTracks
    //getTrackId(0)
    //getTrackId(1)

    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;

}

