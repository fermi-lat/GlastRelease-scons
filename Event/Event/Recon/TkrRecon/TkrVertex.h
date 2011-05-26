/** @file TkrVertex.h
* @author The Tracking Software Group
*
* $Header$
*/

#ifndef __TkrVertex_H
#define __TkrVertex_H 1

#include <vector>
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/IInterface.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"


static const CLID& CLID_TkrVertex = InterfaceID("TkrVertex", 1, 1);

namespace Event { //Namespace

/** 
* @class TkrVertex
*
* @brief TDS Data Object containing the information for a single tracker reconstructed vertex
*
* This class defines the basic TkrRecon structure containing output information
* for a vertex reconstructed in the LAT Tracker from TkrTracks. This class 
* has been modelled after TkrTrack. 
* Information included is:
* 1) Reconstruction Status information (see StatusBits)
* 2) Summary reconstruction information including:
*    a) Vertex ordering information
*    b) estimated energy, ...
*    c) Variables to help estimate vertex quality
*    d) etc.
* 3) Also created is a TkrVertexTab to establish the relationship (ownership)
*    of tracks by vertices and visa-versa
*
*/

    class TkrVertex : virtual public ContainedObject 
    {   
    public:
        TkrVertex(): m_statusBits(0), m_energy(0.),  m_position(), m_direction(),
                     m_chiSquare(0.), m_quality(0.), m_arcLen1(0.), m_arcLen2(0.),
                     m_doca(0.), m_radlen(0.), m_vtxID(0,0,0,false)  {}

        TkrVertex(idents::TkrId tkrID, double energy, double quality, double chisq, 
                   double rad_len, double doca, double s1, double s2, double z,
                   TkrTrackParams params);
        virtual ~TkrVertex() {}

        //! Retrieve pointer to class defininition structure
        virtual const CLID& clID() const   { return TkrVertex::classID(); }
        static const CLID& classID()       { return CLID_TkrVertex; }

        // Status word bits organized like:
        //        |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
        //                     Ghost                   [Track Topology ] [Track composition]
        enum StatusBits {
            ONETKRVTX         = 0x0001,  //Set if single track vertex
            TWOTKRVTX         = 0x0002,  //Set if 2 track vertex
            WIDEFIRSTCLUSTER  = 0x0004,  //Set if the first cluster is "too wide"
            NEUTRALVTX        = 0x0008,  //Set if vertex includes neutral energy vector
            DOCAVTX           = 0x0010,  //Set if vertex location set by DOCA point
            FIRSTHIT          = 0x0020,  //Set if two tracks share first hit
            STAGVTX           = 0x0040,  //Set if tracks don't start in same plane (staggered)
            CROSSTKR          = 0x0080,  //Set if DOCA location lies inside track hits

            GHOST             = 0x1000   //Set if at least one track in the vertex is a ghost
        };

        /// Utility 
        std::ostream& fillStream( std::ostream& s ) const;

        /// Access to primary quantities on track quality and scattering info
        inline unsigned int getStatusBits()          const {return m_statusBits;}
        inline double       getChiSquare()           const {return m_chiSquare;}
        inline double       getQuality()             const {return m_quality;}


        /// Access to fit specific information
        inline Point        getPosition()            const {return m_position;}
        inline Vector       getDirection()           const {return m_direction;}
        inline double       getEnergy()              const {return m_energy;}
        TkrTrackParams      getVertexParams()        const {return m_params;} 

        /// Access to other Geometry Information
        inline const idents::TkrId getTkrId()        const {return m_vtxID;}
        inline double       getAddedRadLen()         const {return m_radlen;}
        inline double       getTkr1ArcLen()          const {return m_arcLen1;}
        inline double       getTkr2ArcLen()          const {return m_arcLen2;}
        inline double       getDOCA()                const {return m_doca;}

        /// Set functions for initializing the data members
        inline void   setPosition(const Point x)          {m_position   = x;}
        inline void   setDirection(const Vector x)        {m_direction  = x;}
        inline void   setEnergy(double x)                 {m_energy     = x;}
        inline void   setChiSquare(double x)              {m_chiSquare  = x;}
        inline void   setQuality(double x)                {m_quality    = x;}
        inline void   setAddedRadLen(double x)            {m_radlen     = x;}
        inline void   setTkr1ArcLen(double x)             {m_arcLen1    = x;}
        inline void   setTkr2ArcLen(double x)             {m_arcLen2    = x;}
        inline void   setDOCA(double x)                   {m_doca       = x;}
        inline void   setTkrID(idents::TkrId& tkrID)      {m_vtxID      = tkrID;}
        inline void   setParams(TkrTrackParams& params)   {m_params     = params;}
        inline void   setStatusBit(unsigned int status)   {m_statusBits |= status;}
        inline void   clearStatusBits(unsigned int bits= 0xffffffff)            
                                                          {m_statusBits &= ~bits;}

        // Add tracks to the list
        void addTrack(const Event::TkrTrack* pTrack)      {m_tracks.push_back(pTrack);}

        // delete last track in the list
        void deleteTrack() {m_tracks.pop_back();}

        // How many tracks in the vertex?
        int  getNumTracks() const {return m_tracks.size();}

        // Pointers to track info
        SmartRefVector<TkrTrack>::const_iterator getTrackIterBegin() const {return m_tracks.begin();}
        SmartRefVector<TkrTrack>::const_iterator getTrackIterEnd()   const {return m_tracks.end();}

        /// Utilities 
        void writeOut(MsgStream& log) const; 

    private:    
        /// Status
        unsigned int m_statusBits;

        /// The energy, position and direction
        double       m_energy;              // energy associated with vertex
        Point        m_position;            // position of vertex
        Vector       m_direction;           // direction of gamma causing pair conversion vertex

        // Vertex quality information
        double       m_chiSquare;           // Spacial chi-square for combining tracks
        double       m_quality;             // Vertex "Quality" derived from topology & chisq.
        double       m_arcLen1;             // Signed distance from head of track 1 to VTX
        double       m_arcLen2;             // Signed distance from head of track 1 to VTX
        double       m_doca;                // Distance between tracks at VTX location
        double       m_radlen;              // Integrated radiation lengths from end of track 1
                                            // to vertex location
        idents::TkrId  m_vtxID;             // Complete TkrId identifying the details of this vertex
                                            // (This is the TkrId of the first hit after the vertex)
        TkrTrackParams m_params;            // Parameter structure for vertex (includes cov. matrix)
   
        SmartRefVector<TkrTrack> m_tracks;  // List of tracks associated with vertex
    };

    //typedef for the Container
    typedef ObjectVector<TkrVertex>      TkrVertexCol;
    typedef TkrVertexCol::const_iterator TkrVertexConPtr;
    typedef TkrVertexCol::iterator       TkrVertexColPtr;

}; //Namespace

#endif
