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

static const CLID& CLID_TkrCluster = InterfaceID("TkrCluster", 1, 0);

namespace Event {

class TkrCluster : virtual public ContainedObject
  {
  public:
    
    /// enumeration of the view of the cluster
    enum view 
      {
        X, /**< cluster measures X */ 
	Y, /**< cluster measures Y */ 
	XY /**< not valid for clusters */
      };
    
  public:
    
      TkrCluster(): m_tkrId(0,0,0,false) {}
    
    /// Constructor with arguments
    /**
     * Construct a cluster with all of its private members set
     * @param tkrId TKR Volume ID of cluster
     * @param istrip0  first strip
     * @param istripf  last strip
     * @param ToT
     */
    TkrCluster(idents::TkrId tkrId,
	       int istrip0, int istripf, Point position,double ToT, int id=-1)
      : m_tkrId(tkrId),m_strip0(istrip0),m_stripf(istripf),
      m_position(position),m_ToT(ToT),m_flag(false), m_id(id)
      {;}
    
    virtual ~TkrCluster() {}
    
    //! Retrieve pointer to class defininition structure
    virtual const CLID& clID() const   { return TkrCluster::classID(); }
    static const CLID& classID()       { return CLID_TkrCluster; }
    
    
    
    // set methods
    /// sets the flag of a cluster
    inline void flag(int flag=1) {m_flag = flag;}
    /// clears the flag of a cluster
    inline void unflag()         {m_flag = 0;}
    
    // get methods
    idents::TkrId getTkrId() const {return m_tkrId;}
    int         tower()      const; //DANGEROUS: SOON TO BE REMOVED
    inline int     id()         const {return m_id;}

    inline int towerX()      const {return m_tkrId.getTowerX();} 
    inline int towerY()      const {return m_tkrId.getTowerY();} 
    inline int  plane()      const {return m_tkrId.getPlane();}
    inline int  layer()      const {return m_tkrId.getLayer();}
    inline view v()          const {return intToView(m_tkrId.getView());}

    inline int  firstStrip() const {return m_strip0;}
    inline int  lastStrip()  const {return m_stripf;}
    inline double ToT()      const {return m_ToT;}
    
    int    chip();
    double strip();
    double size ();
    
    Point position()       const {return m_position;}
    
    /// returns true if the cluster has been flagged
    bool hitFlagged()      const {return (m_flag!=0);}
    
    /// writes out the information of the cluster if msglevel is set to debug
    void writeOut(MsgStream& log) const;
    
    /// converts the view integer to enum view
    static enum view intToView(int);
    /// converts enum to int; guarantees that X->0 and Y->1
    static int viewToInt(view v);
    
  private:
    
    /// volume id
    idents::TkrId m_tkrId;  
    /// initial strip address of the cluster
    int m_strip0;
    /// final strip address of the cluster
    int m_stripf;
    /// space position of the cluster
    Point m_position;
    /// ToT value of the cluster
    double m_ToT;    
    /// flag of the cluster, used during pattern recognition
    int m_flag;

    //TO BE REMOVED
    int m_id;
    
  };

//typedef for the Container (to be stored in the TDS)
//typedef ObjectVector<TkrCluster>       TkrClusterCol;
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
        //return left.getLayer() > right.getLayer();
        return left < right;
    }
};

typedef std::multimap<idents::TkrId,SmartRef<Event::TkrCluster>,CompareTkrIds> clusIdMMap;
typedef std::pair<idents::TkrId, SmartRef<Event::TkrCluster> >                 clusIdPair;

static const CLID& CLID_TkrIdClusterMMap  = InterfaceID("TkrIdClusterMMap",  1, 0);

class TkrIdClusterMMap : public clusIdMMap, public DataObject 
{
    //! Retrieve pointer to class defininition structure
    virtual const CLID& clID() const   { return TkrIdClusterMMap::classID(); }
    static const CLID& classID()       { return CLID_TkrIdClusterMMap; }
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

};

#endif // TKRCLUSTER_H
