/** @file TkrCluster.h
* @author Tracy Usher, Leon Rochester
*
* $Header$

*/
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


#include "GaudiKernel/IInterface.h"

static const CLID& CLID_TkrCluster = InterfaceID("TkrCluster", 7, 0);

namespace Event {

    /** 
    * @class TkrCluster
    *
    * @brief Contains the data members which specify a TKR cluster, 
    * and access methods
    *
    * Adapted from SiCluster of Jose Hernando
    *
    */

    class TkrCluster : virtual public ContainedObject
    {
    public:

        // once we have an official release, the version number can be used
        //  to allow backward compatibility

        //enum {VERSION = 3};

        // fields and shifts of the status word, which together make a mask
        //  

        // Status word is organized as follows:
        //
        // Low-order bits (0-15, right to left):
        //
        // |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
        //
        //            A   A         T   S   G                 R     u  {Cntlr}  U
        //            l   l         o   a   h                 e     s    0=0    s
        //            o   o         T   m   o                 m     e    1=1    e
        //            n   n         2   e   s                 o     d   2=1+2   d
        //            E   e         5   T   t                 v     C      
        //            n             5   k                     e     R   
        //            d                                       d
        //
        // High-order bits (16-31, right to left):
        //
        // |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
        //                                                  
        //        L   P                           C           G                 G
        //        y   l                           o           d                 h
        //        r   n                           m           T                 o
        //        O   O                           p           r                 s
        //        f   f                           o           e                 t
        //        f   f                           s           e                 D
        //                                      [Tree-based bits]               
                                                                            
        enum { 
            fieldUSED        =  1,    // tells whether cluster is used on a track
            fieldEND         =  3,    // identifies controller, 0, 1, 2=mixed
            fieldUSEDCR      =  1,    // cluster used on a CR track
            fieldREMOVED     =  1,    // cluster removed by Ghost Filter
            fieldGHOST       =  1,    // cluster is marked as a ghost
            fieldSAMETRACK   =  1,    // This cluster belongs to a track with a 255 or ghost
            field255         =  1,    // cluster is marked as a ToT==255
            fieldALONE       =  1,    // cluster is alone in its plane
            fieldALONEEND    =  1,    // cluster is alone in its readout end
            fieldDIAGNOSTIC  =  1,    // ghost cluster discovered from TEM diags
            fieldSAMETRACKD  =  1,    // This cluster from a track with a 255 or a diagnostic ghost
            fieldTREEBITS    = 15,    // Reserving this 4 bit field for association to trees
            fieldONAGOODTREE =  1,    // This clusters has been associated to a "good" tree
            fieldCOMPOSITE   =  1,    // This cluster is a "composite" cluster from tree based tracking
            fieldMERGED      =  1,    // This cluster was merged into a super cluster
            filedMERGERESULT =  1,    // This clusters is the result of merging other clusters
            fieldPLANEOFFSET =  1,    // to calculate Plane number from Tray/Face (1 for LAT)
            fieldLAYEROFFSET =  1     // to calculate Layer number from Plane (0 for LAT)
        };
        enum {    
            shiftUSED        =  0,
            shiftEND         =  1,
            shiftUSEDCR      =  3,
            shiftREMOVED     =  4,
            shiftGHOST       =  8,
            shiftSAMETRACK   =  9,
            shift255         = 10,
            shiftALONE       = 12,
            shiftALONEEND    = 13,
            shiftDIAGNOSTIC  = 16,
            shiftSAMETRACKD  = 17,
            shiftTREEBITS    = 20,
            shiftCOMPOSITE   = 23,
            shiftMERGED      = 24,
            shiftMERGERESULT = 25,
            shiftPLANEOFFSET = 29,
            shiftLAYEROFFSET = 30 
        };
        enum maskType {
            maskUSED        = fieldUSED<<shiftUSED,
            maskEND         = fieldEND<<shiftEND,
            maskUSEDCR      = fieldUSEDCR<<shiftUSEDCR,
            maskREMOVED     = fieldREMOVED<<shiftREMOVED,
            maskGHOST       = fieldGHOST<<shiftGHOST,
            maskSAMETRACK   = fieldSAMETRACK<<shiftSAMETRACK,
            mask255         = field255<<shift255,
            maskALONE       = fieldALONE<<shiftALONE,
            maskALONEEND    = fieldALONEEND<<shiftALONEEND,
            maskDIAGNOSTIC  = fieldDIAGNOSTIC<<shiftDIAGNOSTIC,
            maskSAMETRACKD  = fieldSAMETRACKD<<shiftSAMETRACKD,
            maskTREEBITS    = fieldTREEBITS<<shiftTREEBITS,
            maskONAGOODTREE = fieldONAGOODTREE<<shiftTREEBITS,
            maskCOMPOSITE   = fieldCOMPOSITE<<shiftCOMPOSITE,
            maskMERGED      = fieldMERGED<<shiftMERGED,
            maskMERGERESULT = filedMERGERESULT<<shiftMERGERESULT,
            maskPLANEOFFSET = fieldPLANEOFFSET<<shiftPLANEOFFSET,
            maskLAYEROFFSET = fieldLAYEROFFSET<<shiftLAYEROFFSET,
            maskZAPGHOSTS   = mask255|maskGHOST|maskSAMETRACK|maskDIAGNOSTIC|maskSAMETRACKD,
            maskUSEDANY     = maskUSED|maskUSEDCR
        };

        TkrCluster(): m_tkrId(0,0,0,false) {}

        /// Constructor with arguments
        /**
        * Construct a cluster with all of its private members set
        * @param tkrId TKR Volume ID of cluster
        * @param istrip0  first strip
        * @param istripf  last strip
        * @param nBad number of bad strips in the cluster
        * @param position global position of cluster center
        * @param rawToT raw ToT from TkrDigi
        * @param ToT corrected ToT
        * @param status word containting various useful info
        */
        TkrCluster(idents::TkrId tkrId, int istrip0, int istripf, 
            Point position, int rawToT, float ToT, unsigned int status, int nBad)
            : m_tkrId(tkrId),m_strip0(istrip0),m_stripf(istripf), m_nBad(nBad),
            m_position(position), m_rawToT(rawToT), m_ToT(ToT), m_status(status) 

        { }

        virtual ~TkrCluster() {}

        //! Retrieve pointer to class defininition structure
        virtual const CLID& clID() const   { return TkrCluster::classID(); }
        static const CLID& classID()       { return CLID_TkrCluster; }

        /// set the corrected ToT
        void setMips(float ToT) {m_ToT = ToT;}

        // get methods
        idents::TkrId getTkrId()   const {return m_tkrId;}
        int           tower()      const {return idents::TowerId(m_tkrId.getTowerX(),
            m_tkrId.getTowerY()).id();} //DANGEROUS: SOON TO BE REMOVED??

        /// sets the USED flag of a cluster (or clear if flag==0), legacy
        inline void flag(int flag=1) { setStatusBits(maskUSED, flag);}
        /// clears the used flag of a cluster, legacy
        inline void unflag() { clearStatusBits(maskUSED); }

        inline int    firstStrip() const {return m_strip0;}
        inline int    lastStrip()  const {return m_stripf;}
        inline double strip()      const {return 0.5*(m_strip0+m_stripf);}
        inline double size ()      const {return std::abs(m_stripf-m_strip0) + 1.;}

        inline double ToT()        const {return getRawToT();}

        inline Point position()    const {return m_position;}
        inline int   getStatusWord() const {return m_status;}
        inline int   getNBad()     const {return m_nBad;}

        /// construct plane from tray/face
        inline int getPlane() const {
            return m_tkrId.getTray()*2 + m_tkrId.getBotTop() - 1 + getPlaneOffset(); 
        }
        // construct layer from Plane
        inline int getLayer() const { 
            return (getPlane() + getLayerOffset())/2 ; }
        // cluster used on a standard track, legacy
        inline bool hitFlagged()     const { return isSet(maskUSED); }
        // returns chip number, hardwired strips/chip = 64
        inline int    chip()  const { return m_strip0/64;}

        /// writes out the information of the cluster if msglevel is set to debug
        inline std::ostream& fillStream( std::ostream& s ) const;

        /// retrieves raw ToT (will be raw or corrected depending on the version)
        inline double getRawToT() const { 
            return ( m_rawToT);
        }
        /// retrieve corrected ToT (zero for old-style records)
        inline double getMips() const {
            return  m_ToT ;
        }
        /// retrieve end
        inline int getEnd() const {
            return ((m_status&maskEND)>>shiftEND);
        }

        /// set/clear/test any bits (including USED)
        inline bool isSet(unsigned int mask) const { return (m_status&mask)>0; }
        inline void setStatusBits(unsigned int mask, int flag=1) {
            if(flag!=0) { m_status |= mask;  }
            else        { clearStatusBits(mask); }
        }
        inline void clearStatusBits(unsigned int mask) {
            m_status &= ~mask;
        }

        /// set USEDCR bit
        inline void setUSEDCRBit() {
            m_status |= maskUSEDCR;
        }

    private:

        inline int getPlaneOffset() const { 
            return  ((m_status&maskPLANEOFFSET)>>shiftPLANEOFFSET);}
        inline int getLayerOffset() const { 
            return ((m_status&maskLAYEROFFSET)>>shiftLAYEROFFSET); }

        /// tracker volume id
        idents::TkrId m_tkrId;  
        /// initial strip address of the cluster
        int    m_strip0;
        /// final strip address of the cluster
        int    m_stripf;
        /// number of bad strips in this cluster
        int m_nBad;
        /// space position of the cluster
        Point  m_position;
        /// raw ToT
        int m_rawToT;
        /// Corrected ToT value of the cluster
        float m_ToT;    
        /// various odds and ends
        unsigned int   m_status;
    };

    std::ostream& TkrCluster::fillStream( std::ostream& s ) const
    {
        // Purpose: writes out debug info
        // Inputs:  message stream
        // Outputs: data written to message stream

        s << " tray " << getTkrId().getTray() 
            << " face " << getTkrId().getBotTop()
            << " XY " << getTkrId().getView()
            << " pos (" << m_position.x() << ", " << m_position.y()
            << ", " << m_position.z() << ") "
            << " i0-if " << m_strip0 <<"-"<< m_stripf
            << " rawToT " << m_rawToT
            << " Mips " << m_ToT;
        //<< " status " << m_status; // std::hex << m_status << std::dec;
        return s;
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
