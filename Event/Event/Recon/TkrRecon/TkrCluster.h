//  associated with Tracker for all the evt Status
//
//---------------------------------------------------

#ifndef TKRCLUSTER_H
#define TKRCLUSTER_H 1

//#include <vector>
#include <map>
#include <set>
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/MsgStream.h"

#include "idents/TkrId.h"
#include "idents/TowerId.h"
#include "geometry/Point.h"

/** 
* @class TkrCluster
*
* @brief Contains the data members which specify a TKR cluster, 
* and access methods
*
* Adapted from SiCluster of Jose Hernando
*
* @author Tracy Usher, Leon Rochester
*
* $Header$
*/

#include "GaudiKernel/IInterface.h"

static const CLID& CLID_TkrCluster = InterfaceID("TkrCluster", 2, 0);

namespace Event {

    class TkrCluster : virtual public ContainedObject
    {
    public:

        // once we have an official release, the version number can be used
        //  to allow backward compatibility

        enum {VERSION = 2};

        // fields and shifts of the status word, which together make a mask
        //  

        enum { 
            fieldUSED        = 1,    // tells whether cluster is used on a track
            fieldEND         = 3,    // identifies controller, 0, 1, 2=mixed
            fieldRAWTOT      = 0xff, // raw ToT
            fieldPLANEOFFSET = 1,    // to calculate Plane number from Tray/Face (1 for LAT)
            fieldLAYEROFFSET = 1,    // to calculate Layer number from Plane (0 for LAT)
            fieldVERSION     = 0xf   // version of class, room for 15!
        };
        enum {    
            shiftUSED        =  0,
            shiftEND         =  1,
            shiftRAWTOT      =  8,
            shiftPLANEOFFSET = 26,
            shiftLAYEROFFSET = 27,
            shiftVERSION     = 28 
        };
        enum {
            maskUSED        = fieldUSED<<shiftUSED,
            maskEND         = fieldEND<<shiftEND,
            maskRAWTOT      = fieldRAWTOT<<shiftRAWTOT,
            maskPLANEOFFSET = fieldPLANEOFFSET<<shiftPLANEOFFSET,
            maskLAYEROFFSET = fieldLAYEROFFSET<<shiftLAYEROFFSET,
            maskVERSION     = fieldVERSION<<shiftVERSION
        };

        TkrCluster(): m_tkrId(0,0,0,false) {}

        /// Constructor with arguments
        /**
        * Construct a cluster with all of its private members set
        * @param tkrId TKR Volume ID of cluster
        * @param istrip0  first strip
        * @param istripf  last strip
        * @param position global position of cluster center
        * @param ToT corrected ToT
        */
        TkrCluster(idents::TkrId tkrId, int istrip0, int istripf, Point position, 
            double ToT, unsigned int status)
            : m_tkrId(tkrId),m_strip0(istrip0),m_stripf(istripf),
            m_position(position), m_ToT(ToT), m_status(status), m_id(-1)
        {}

        virtual ~TkrCluster() {}

        //! Retrieve pointer to class defininition structure
        virtual const CLID& clID() const   { return TkrCluster::classID(); }
        static const CLID& classID()       { return CLID_TkrCluster; }

        // get methods
        idents::TkrId getTkrId()   const {return m_tkrId;}
        int           tower()      const {return idents::TowerId(m_tkrId.getTowerX(),
            m_tkrId.getTowerY()).id();} //DANGEROUS: SOON TO BE REMOVED??

        inline int    firstStrip() const {return m_strip0;}
        inline int    lastStrip()  const {return m_stripf;}
        inline double strip()      const {return 0.5*(m_strip0+m_stripf);}
        inline double size ()      const {return std::abs(m_stripf-m_strip0) + 1.;}

        inline double ToT()        const {return m_ToT;}

        inline Point position()    const {return m_position;}
        inline int   getStatusWord() const {return m_status;}

        /// construct plane from tray/face
        inline int getPlane() const {
            return m_tkrId.getTray()*2 + m_tkrId.getBotTop() - getPlaneOffset(); 
        }
        // construct layer from Plane
        inline int getLayer() const { 
            return (getPlane() + getLayerOffset())/2 ; }
        // cluster used on a track
        bool hitFlagged()     const { 
            return ((maskVERSION&m_status)!=0 ? ((m_status&maskUSED)>0) : m_status!=0);}
        // returns chip number, hardwired strips/chip = 64
        inline int    chip()  const { return m_strip0/64;}

        /// writes out the information of the cluster if msglevel is set to debug
        inline void writeOut(MsgStream& log) const;

        // set methods
        /// sets the used flag of a cluster
        inline void flag(int flag=1) {
            unflag();
            if(flag>0) m_status |= maskUSED;
        }
        /// clears the used flag of a cluster
        inline void unflag()         {m_status &= ~maskUSED;}

        /// Everything below here is a candidate for removal
        //inline int    id()         const {return m_id;}

    private:

        inline int getPlaneOffset() const { 
            return ((maskVERSION & m_status)!=0 ? ((m_status&maskPLANEOFFSET)>>shiftPLANEOFFSET) : 1); }
        inline int getLayerOffset() const { 
            return ((maskVERSION & m_status)!=0 ?((m_status&maskLAYEROFFSET)>>shiftLAYEROFFSET) : 0); }

        /// tracker volume id
        idents::TkrId m_tkrId;  
        /// initial strip address of the cluster
        int    m_strip0;
        /// final strip address of the cluster
        int    m_stripf;
        /// space position of the cluster
        Point  m_position;
        /// ToT value of the cluster
        float m_ToT;    
        /// various odds and ends
        unsigned int   m_status;

        /// Everything  below here is a candidate for removal
        //TO BE REMOVED
        int m_id;
    };

    void TkrCluster::writeOut(MsgStream& log) const
    {
        // Purpose: writes out debug info
        // Inputs:  message stream
        // Outputs: data written to message stream

        log << MSG::DEBUG;
        if (log.isActive()) {
            log << " tray " << getTkrId().getTray() 
                << " face " << getTkrId().getBotTop()
                << " XY " << getTkrId().getView()
                << " xpos  " << m_position.x()  << " ypos   " << m_position.y()
                << " zpos  " << m_position.z()
                << " i0-if " << m_strip0 <<"-"<< m_stripf;
        }
        log << endreq;
    }


    //typedef for the Container (to be stored in the TDS)
    typedef ObjectVector<TkrCluster>       TkrClusterCol;
    typedef TkrClusterCol::const_iterator  TkrClusterColConItr;
    typedef TkrClusterCol::iterator        TkrClusterColItr;

    // typedef for creating a vector in the main track object
    typedef SmartRefVector<TkrCluster>     TkrClusterVec;
    typedef TkrClusterVec::const_iterator  TkrClusterVecConItr;
    typedef TkrClusterVec::iterator        TkrClusterVecItr;

    // typedefs for mapping TkrIds to clusters
    struct CompareTkrIds
    {
    public:
        bool operator()(const idents::TkrId left, const idents::TkrId right) const
        {
            return left < right;
        }
    };

    typedef std::map<idents::TkrId,Event::TkrClusterVec,CompareTkrIds> clusVecIdMap;
    typedef std::pair<idents::TkrId, Event::TkrClusterVec >            clusVecIdPair;

    static const CLID& CLID_TkrIdClusterMap  = InterfaceID("TkrIdClusterMap",  1, 0);

    class TkrIdClusterMap : public clusVecIdMap, public DataObject 
    {
        //! Retrieve pointer to class defininition structure
        virtual const CLID& clID() const   { return TkrIdClusterMap::classID(); }
        static const CLID& classID()       { return CLID_TkrIdClusterMap; }
    };
}

#endif // TKRCLUSTER_H
