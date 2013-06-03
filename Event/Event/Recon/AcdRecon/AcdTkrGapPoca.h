#ifndef EVENT_ACDTKRGAPPOCA_H
#define EVENT_ACDTKRGAPPOCA_H

#include "AcdTkrLocalData.h"
#include "AcdPocaData.h"

#include "idents/AcdGapId.h"
#include <vector>

#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"

static const CLID& CLID_AcdTkrGapPocaCol = InterfaceID("AcdTkrGapPocaCol", 1, 0);

/** 
 * @class Event::AcdTkrGapPoca
 * @brief TDS object which stores information about the Point of Closest Approach (POCA) between 
 * an extrapolated track and a Gap in the ACD
 *  
 * Most of the structure of the object comes from the base classes AcdTkrLocalCoords and AcdPocaData
 * 
 * The class adds only enough information to define the involved in the POCA.
 *    - const idents::AcdGapId& getId()  
 *      - which returns the ID of the gap
 *    - const int getTrackIndex()  
 *      - which returns the index of the track which did the hitting
 *
 * 
 * @author Eric Charles
 *
 * $Header$
 */

namespace Event {

  class AcdTkrGapPoca : public Event::AcdTkrLocalCoords, public Event::AcdPocaData  {

  typedef HepGeom::Point3D<double> HepPoint3D;
  
  private:

   static unsigned s_uid;
    
  public:
    
    /// Default constructor.  
    AcdTkrGapPoca();
    
    /// Constructor for use in transient -> persistent conversion 
    /// Takes arguements as they are stored in ROOT
    AcdTkrGapPoca(const idents::AcdGapId& acdId, int trackIndex, 
                  const float active[2], 
                  float vetoSigmaHit, float vetoSigmaProj, float vetoSigmaProp,
                  int volumePlane, float arcLengthToPlane, float cosTheta, 
                  const HepPoint3D& global, const float localPosition[2], 
                  const CLHEP::HepSymMatrix& localCovProj, const CLHEP::HepSymMatrix& localCovProp,
                  int volume, int region, float arcLength, 
                  float doca, float docaErrProj, float docaErrProp,
                  const Point& poca, const Vector& voca);

    /// Old Constructor for backwards compatiblity
    AcdTkrGapPoca( const idents::AcdGapId& acdId, int trackIndex, 
                   const Event::AcdTkrLocalCoords& local, const Event::AcdPocaData& pocaData );

  
    /// Copy constructor
    AcdTkrGapPoca(const Event::AcdTkrGapPoca& params);
    
    /// Assignment operator
    Event::AcdTkrGapPoca& operator=(const Event::AcdTkrGapPoca& params);

    /// Equality test operator, requires identity
    bool operator==(const Event::AcdTkrGapPoca& other) const {
      return this == &other;
    }

    /// Comparison operator, sorts by vetoSigma
    bool operator<(const Event::AcdTkrGapPoca& other) const;     

    /// Return the AcdId of the Gap 
    inline const idents::AcdGapId& getId() const { return m_id; }
    
    /// Return the index of the associated track
    inline int trackIndex() const {
      return m_trackIndex;
    }

    /// combine the sigma from the hit with the sigma from the track
    inline float vetoSigma2() const {
      return ((m_vetoSigmaHit*m_vetoSigmaHit) + (m_vetoSigmaProj*m_vetoSigmaProj));
    }
    
    /// An estimator of the the chance that a track in this gap would produce no signal
    inline float vetoSigmaHit() const {
      return m_vetoSigmaHit;
    }

    /// An estimator of the chance that a track went into this gap
    inline float vetoSigmaProj() const {
      return m_vetoSigmaProj;
    }

    /// An estimator of the chance that a track went into this gap
    inline float vetoSigmaProp() const {
      return m_vetoSigmaProp;
    }

    /// return the UID for sorting
    inline unsigned getUID() const {
      return m_uid;
    }

    /// set only the data at this level
    inline void set(const idents::AcdGapId& gapId, int trackIndex, 
                    float vetoSigmaHit, float vetoSigmaProj, float vetoSigmaProp) {
      m_id = gapId; m_trackIndex = trackIndex; 
      m_vetoSigmaHit = vetoSigmaHit; m_vetoSigmaProj = vetoSigmaProj; m_vetoSigmaProp = vetoSigmaProp;
    }

    /// reset all the values to their default
    virtual void ini();
    
    /// Print out this structure
    virtual void writeOut(MsgStream& stream) const;
    
  private:
  
    /// The ID of the Gap tile
    idents::AcdGapId m_id;
    
    /// The index of the associated track
    int m_trackIndex;
    
    /// An estimator of the chance that a track in this gap would produce no signal
    float m_vetoSigmaHit;

    /// An estimator of the chance that a track went into this gap
    float m_vetoSigmaProj;

    /// An estimator of the chance that a track went into this gap
    float m_vetoSigmaProp;

     /// A UID for sorting
     unsigned m_uid;
  };
 

  /*! 
   * @class AcdTkrGapPocaCol
   *
   *  @brief TDS class a collection of AcdTkrGapPoca objects
   *  
   * It inherits from DataObject
   * (to be a part of Gaudi TDS) and from std::vector of pointers
   * to AcdTkrGapPoca objects. Some methods just rename the std::vector
   * methods with the same functionality for backward compartibility.
   *
   *
   * @author Eric Charles
   *
   * @todo replace this class by typedef to ObjectVector, make corresponding
   *       changes in AcdReconAlg
   */
    
  class AcdTkrGapPocaCol : public DataObject, public std::vector<AcdTkrGapPoca*> {
    
  public:
    
    AcdTkrGapPocaCol(const std::vector<AcdTkrGapPoca*>& acdTkrIntersections);

    AcdTkrGapPocaCol() { clear();}
    
    /// destructor - deleting the clusters pointed
    /// by the vector elements
    ~AcdTkrGapPocaCol() { del();}
        
    
    // GAUDI members to be use by the converters
    static const CLID& classID() {return CLID_AcdTkrGapPocaCol;}
    virtual const CLID& clID() const {return classID();}
    
    /// add new cluster
    void add(AcdTkrGapPoca* cl) {push_back(cl);}
    
    /// get the number of clusters in collection
    int num()                  const {return size();}
    
    /// get pointer to the cluster with given number 
    AcdTkrGapPoca* get(int i) const {return operator[](i);}
    
    /// delete all intersections pointed by the vector elements
    void del();
    
    /// write information for all clusters to the ascii file 
    /// for debugging purposes
    virtual void writeOut(MsgStream& stream) const;
    
  protected:
    
    /// does the same function as clear() 
    virtual void ini();
        
  };

  

};
#endif
