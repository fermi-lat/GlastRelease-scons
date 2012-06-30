#ifndef CalCluster_H
#define CalCluster_H

#include <vector>
#include <string>
#include <map>
#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/MsgStream.h"
#include "Event/RelTable/RelTable.h"
#include "Event/Recon/CalRecon/CalXtalsParams.h"
#include "Event/Recon/CalRecon/CalFitParams.h"
#include "Event/Recon/CalRecon/CalParams.h"
#include "Event/Recon/CalRecon/CalMomParams.h"
#include "Event/Recon/CalRecon/CalClassParams.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalMSTreeParams.h"

#define NUMCALLAYERS 8

static const CLID& CLID_CalClusterCol = InterfaceID("CalClusterCol", 1, 0);

/**
*  @class CalCluster
*
*
*  @brief Defines a CalRecon TDS object which will contain the results of
*         clustering - the association of CalXtalRecData objects. This includes
*         defining a RelTable relating clusters with the CalXtalRecData
*         objects comprising them. 
* 
*         Main changes from old CalCluster.h:
*         1) Remove all "final" energy recon information (leakage, profile fitting, etc.)
*            This will be moved to the higher level object - CalEventEnergy
*         2) Incorporate layer data into a subclass and store in a vector. The idea is 
*            that the vector will always have 8 "empty" elements, presuming that external
*            usage is to loop over Cal "layers".
*         3) Provide individal "set" and "get" routines (get was there, now include set)
*  
*  @author The CalRecon Rewrite Group
*
* $Header$
*/

namespace Event { //Namespace Event

  /**
   *   Define a subclass which will hold data by later
   *
   */
  class CalClusterLayerData
  {
  public:
    CalClusterLayerData() : m_energy(0.), m_position(0.,0.,0.), m_rmsSpread(0.,0.,0.) {}
      CalClusterLayerData(double energy, Point position, Vector rmsSpread) :
        m_energy(energy), m_position(position), m_rmsSpread(rmsSpread) {}
        
    ~CalClusterLayerData() {}

    double         getEnergy()    const {return m_energy;}
    const  Point&  getPosition()  const {return m_position;}
    const  Vector& getRmsSpread() const {return m_rmsSpread;}
    
  private:
    
    /// Energy deposition in crystals in this cluster and layer.
    double m_energy;
    // Average position in this layer.
    Point  m_position;
    // Quadratic position spread for this layer.
    Vector m_rmsSpread;
  };
  
  typedef std::vector<CalClusterLayerData> CalClusterLayerDataVec;
  
  
  /// Note: when changing this class, the following files should be changed accordingly:
  /// - RootConvert/src/Recon/CalClusterConvert.cxx
  ///   
  /// - reconRootData/reconRootData/CalCluster.h
  /// - reconRootData/src/CalCluster.cxx
  ///   all the interfaces and members must be compatible with the TDS version.
  /// 
  
  class CalCluster : public CalClusterLayerDataVec, virtual public ContainedObject 
  {     
  public:
  
    /// Status word bits organized like:
    /// low:   |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
    ///         < Not used   >    <Analysis Status> <   Found By   >  < Cluster Type >
    /// high:  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
    ///         <  Not yet used                                                        >
    enum StatusBits {
      /// @param ZERO
      ZERO            = 0x00000000, 
      /// @param ALLXTALS Bit set means all xtals in the event are here - super set
      ALLXTALS        = 0x00000001, 
      // DC: redundant with ALLXTALS
      //    /// @param ISOLATEDCLUSTER Bit set means this is a subset of all Xtals
      //                     ISOLATEDCLUSTER = 0x00000002,  
      // DC: not used ?
      //        /// @param MIPTRACK Bit set means this is a subset of all Xtals associate with MIP
      //                     MIPTRACK        = 0x00000004,  
      // DC: to be replaced with a string
      //        /// @param SIMPLECLUSTER Bit set means cluster found using SimpleClusterTool
      //                     SIMPLECLUSTER   = 0x00000010,  
      //        /// @param FUZZYCLUSTER Bit set means cluster found using FuzzyClusterTool
      //                     FUZZYCLUSTER    = 0x00000020,
      //        /// @param MIPCLUSTER Bit set means cluster found using CalMIPFinderTool
      //                     MIPCLUSTER      = 0x00000040,
      // DC: CENTROID SHOULD BE REPLACED WITH AXIS
      /// @param CENTROID Bit set means centroid found
      CENTROID        = 0x00000100, 
      /// @param MOMENTS Bit set means moments analysis run
      MOMENTS         = 0x00000200,
      // DC: not used ?
      //        /// @param ENERGYCORR Bit set means one or more energy corrections were run
      //                         ENERGYCORR      = 0x00000400,
      // DC: not used ?
      //        /// @param MIPFIT Bit set means a MIP like track was fitted
      //                         MIPFIT          = 0x00000800
      /// @param MSTTREE Bit set means a Minimum Spanning Tree information is available
      MSTTREE          = 0x0001000,
      /// @param CLASSIFIED Bit set means clusters have been classified
      CLASSIFIED       = 0x0002000
      
    };
    
    /// No parameters constructor.
    CalCluster(int numCalLayers = NUMCALLAYERS) :
      CalClusterLayerDataVec(numCalLayers)
    { iniCluster(); }

    /// It's not enrirely clear to me why we don't have a constructor setting
    /// all the class members, but I'll leave things as they are for the time
    /// being (Luca Baldini, Dec 22, 2010).
       
    /// Destructor.
    virtual ~CalCluster() {}
    
    /// Initialize an existing object from all parameters.
    void initialize(const CalXtalsParams& xtalsParams,
                    const CalMSTreeParams& mstParams,
                    const CalFitParams& fitParams,
                    const CalMomParams& momParams,
                    const CalClassParams& classParams);

    /// Access methods to the main objects.
    inline const std::string & getProducerName()     const { return m_producerName; }
    inline unsigned int getStatusBits()              const { return m_statusBits; }
    const CalXtalsParams& getXtalsParams()           const { return m_xtalsParams; }
    const CalMSTreeParams& getMSTreeParams()         const { return m_mstParams; }
    const CalFitParams& getFitParams()               const { return m_fitParams; }
    const CalMomParams& getMomParams()               const { return m_momParams; }
    const CalClassParams& getClassParams()           const { return m_classParams; }

    /// Set methods.
    CalMomParams& getMomParamsRef() { return m_momParams; }
    inline void setProducerName(const std::string & producerName)
      { m_producerName = producerName ; }
    inline void setStatusBits( unsigned int statusBits )   { m_statusBits = statusBits; }
    void setXtalsParams(const CalXtalsParams& xtalsParams) { m_xtalsParams = xtalsParams; }
    void setMSTreeParams(const CalMSTreeParams& mstParams) { m_mstParams = mstParams; }
    void setFitParams(const CalFitParams& fitParams)       { m_fitParams = fitParams; }
    void setMomParams(const CalMomParams& momParams)       { m_momParams = momParams; }
    void setClassParams(const CalClassParams& classParams) { m_classParams = classParams; }

    /// Manipulate status bits.
    inline void setStatusBit( StatusBits bitToSet )        { m_statusBits |=  bitToSet ; }
    inline void clearStatusBit( StatusBits bitToClear )    { m_statusBits &= ~bitToClear ; }
    inline bool checkStatusBit( StatusBits bitToCheck ) const
      { return ((m_statusBits&bitToCheck)!=ZERO) ; }

    // A few convenience methods.
    int getNumXtals()          const { return m_xtalsParams.getNumXtals(); }
    int getNumSaturatedXtals() const { return m_xtalsParams.getNumSaturatedXtals(); }
    int getNumTruncXtals()     const { return m_xtalsParams.getNumTruncXtals(); }

    /// Std output facility.
    std::ostream& fillStream(std::ostream& s) const;
    friend std::ostream& operator<< (std::ostream& s, const CalCluster& obj)
    {
      return obj.fillStream(s);
    }

    
    /// The rest of this stuff is mostly obsolete and we should probably get rid of
    /// it, altough I am not sure wheter it would screw somebody else's life up.
    /// I'll leave the cleanup for the next round (Luca Baldini, Dec 22 2010).

    /// This is not to break anything with the new container CalMomParams
    /// It should probably print out something ("obsolete?").
    const CalMomParams& getParams()              const { return m_momParams; }

    // TBD Change names according to the CalMomParams class? Maybe just get rid of them all.
    double getRmsLong()                                 const { return m_momParams.getLongRms(); }
    double getRmsLongAsym()                      const { return m_momParams.getLongRmsAsym(); } 
    double getRmsTrans()                         const { return m_momParams.getTransRms(); }
    double getSkewnessLong()                     const { return m_momParams.getLongSkewness(); }

    // TBD make obsolete in favour of getCentroid()? Or get rid of it, as at this
    // point the code calling this should have the responsibiliy to choose
    // explicitely the centroid it's interested in (i.e. from the fit or from the
    // moments analysis).
    const Point & getPosition()                  const { return m_momParams.getCentroid(); }
    Point getCorPosition(Vector vaxis)     { return m_momParams.getCorCentroid(vaxis); }

    // TBD make obsolete in favour of getAxis()? See the comments three lines above.
    const Vector & getDirection()                const { return m_momParams.getAxis(); }

    /// Write some of CalCluster data to the ASCII output file for debugging purposes
    /// Is this really needed? Or the overload of << is enough?
    void writeOut(MsgStream& stream) const;

        
private:

    /// Reset CalCluster data.
    void iniCluster();
    /// Name of the producer.
    std::string m_producerName;
    /// Status Bits.
    unsigned int m_statusBits;
    /// Generic properties of the xtal collection.
    CalXtalsParams m_xtalsParams;
    /// Output of the Minimum Spanning Tree clustring algorithm.
    CalMSTreeParams m_mstParams;
    /// Output of the fit to the cluster centroid/direction.
    CalFitParams m_fitParams;
    /// Output of the moment analysis.
    CalMomParams m_momParams;
    /// Output of the cluster classification 
    CalClassParams m_classParams;
};


  // typedef for the Gaudi TDS Container.
  // This is the object owning the clusters and responsible for deleting them when
  // it's time to cleanup.
  // We're generally not going to loop over this collection.
  typedef ObjectVector<CalCluster>      CalClusterCol;
  typedef CalClusterCol::iterator       CalClusterColItr;
  typedef CalClusterCol::const_iterator CalClusterColConItr;
  
  // Define the relational table taking us back to CalXtalRecData objects
  typedef Event::RelTable<Event::CalXtalRecData, Event::CalCluster> CalClusterHitTab;
  typedef Event::Relation<Event::CalXtalRecData, Event::CalCluster> CalClusterHitRel;
  typedef RelationList<CalXtalRecData, CalCluster>                  CalClusterHitTabList;
  
}; //Namespace Event

#endif        
