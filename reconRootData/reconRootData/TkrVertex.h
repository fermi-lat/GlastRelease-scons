#ifndef ROOT_TKRVERTEX_H
#define ROOT_TKRVERTEX_H

#include "TObjArray.h"
#include "TVector3.h"
#include "TkrTrackParams.h"
#include "commonRootData/idents/TkrId.h"
#include <vector>
#ifndef R__GLOBALSTL
#ifndef WIN32
using std::vector;
#else
using namespace std;
#endif
#endif


/** 
 * @class TkrVertex
 * @brief Root object representing a reconstructed vertex
 *
 * Includes base summary information about the vertex
 * and an array of track ids associated with the vertex
 * This class should duplicate Event::TkrVertex TDS class
 * 
 * @author Tracy Usher
 *
 * $Header$
 */

class TkrVertex : public TObject
{
public:
    /// Status word bits organized like:
    ///        |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
    ///                                             <Track Topology > <Track composition>
    enum StatusBits {ONETKRVTX = 0x0001,  //Set if single track vertex
                     TWOTKRVTX = 0x0002,  //Set if 2 track vertex
                     MULTKRVTX = 0x0004,  //Set if >2 track vertex
                     DOCAVTX   = 0x0010,  //Set if vertex location set by DOCA point
                     FIRSTHIT  = 0x0020,  //Set if two tracks share first hit
                     STAGVTX   = 0x0040,  //Set if tracks don't start in same plane (staggered)
                     CROSSTKR  = 0x0080}; //Set if DOCA location lies inside track hits
    
    TkrVertex();

    virtual ~TkrVertex();

    void Clear(Option_t *option="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const TkrVertex &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;

    /// Access to primary quantities on track quality and scattering info
    inline UInt_t   getStatusBits()   const {return m_statusBits;}
    inline Double32_t getChiSquare()    const {return m_chiSquare;}
    inline Double32_t getQuality()      const {return m_quality;}


    /// Access to fit specific information
    inline TVector3 getPosition()     const {return m_position;}
    inline TVector3 getDirection()    const {return m_direction;}
    inline Double32_t getEnergy()       const {return m_energy;}
    const TkrTrackParams& getVertexParams() const {return m_params;} 

    /// Access to other Geometry Information
    inline Double32_t   getAddedRadLen()  const {return m_radlen;}
    inline Double32_t   getTkr1ArcLen()   const {return m_arcLen1;}
    inline Double32_t   getTkr2ArcLen()   const {return m_arcLen2;}
    inline Double32_t   getDOCA()         const {return m_doca;}

    inline const commonRootData::TkrId getTkrId() const {return m_vtxID;}

    /// Set functions for initializing the data members
    inline void setPosition(const TVector3 x)          {m_position   = x;}
    inline void setDirection(const TVector3 x)         {m_direction  = x;}
    inline void setEnergy(Double32_t x)                  {m_energy     = x;}
    inline void setChiSquare(Double32_t x)               {m_chiSquare  = x;}
    inline void setQuality(Double32_t x)                 {m_quality    = x;}
    inline void setAddedRadLen(Double32_t x)             {m_radlen     = x;}
    inline void setTkr1ArcLen(Double32_t x)              {m_arcLen1    = x;}
    inline void setTkr2ArcLen(Double32_t x)              {m_arcLen2    = x;}
    inline void setDOCA(Double32_t x)                    {m_doca       = x;}
    inline void setTkrID(commonRootData::TkrId& tkrID) {m_vtxID      = tkrID;}
    inline void setParams(TkrTrackParams& params)      {m_params     = params;}
    inline void setStatusBit(UInt_t status)            {m_statusBits |= status;}
    inline void clearStatusBits()                      {m_statusBits = 0;}

    //! Provide info to any hit on the track
    UInt_t      getNumTracks() const          {return m_tracks.GetEntries(); }
    TObject*    getTrack(UInt_t idx) const    {return m_tracks[idx];   }
    void        addTrack(TObject* track)      {m_tracks.Add(track);}

private:
    
    UInt_t                m_statusBits;    // Status bits to describe vertex

    Double32_t              m_energy;        // energy associated with vertex
    TVector3              m_position;      // position of vertex
    TVector3              m_direction;     // direction of gamma causing pair conversion vertex

    Double32_t              m_chiSquare;     // Spacial chi-square for combining tracks
    Double32_t              m_quality;       // Vertex "Quality" derived from topology & chisq.
    Double32_t              m_arcLen1;       // Signed distance from head of track 1 to VTX
    Double32_t              m_arcLen2;       // Signed distance from head of track 1 to VTX
    Double32_t              m_doca;          // Distance between tracks at VTX location
    Double32_t              m_radlen;        // Integrated radiation lengths from end of track 1
    
    commonRootData::TkrId m_vtxID;         // Complete TkrId identifying the details of this vertex
		                                   // (This is the TkrId of the first hit after the vertex)
    TkrTrackParams        m_params;        // Parameter structure for vertex (includes cov. matrix)

    //! List of track IDs used by this vertex
    TObjArray             m_tracks;

    ClassDef(TkrVertex,3)
};
#endif
