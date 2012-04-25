#ifndef EVENT_ACDTKRHITPOCA_H
#define EVENT_ACDTKRHITPOCA_H

#include "AcdTkrLocalData.h"
#include "AcdPocaData.h"

/** 
 * @class Event::AcdTkrHitPoca
 * @brief TDS object information about the Point of Closest Approach (POCA) between an extrapolated track
 *  and a hit Acd element (tile or ribbon).  
 *  
 * Most of the structure of the object comes from the base classes AcdTkrLocalCoords and AcdPocaData
 * 
 * The class adds only enough information to define the involved in the POCA.
 *    - const idents::AcdId& getId()  
 *      - which returns the ID of the hit element
 *    - const int getTrackIndex()  
 *      - which returns the index of the track which did the hitting
 *
 * @author Eric Charles
 *
 * $Header$
 */

#include <vector>
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"

#include "idents/AcdId.h"

static const CLID& CLID_AcdTkrHitPocaCol = InterfaceID("AcdTkrHitPocaCol", 1, 0);

namespace Event {

  class AcdTkrHitPoca : public AcdTkrLocalCoords, public AcdPocaData  {
    
  public:
    
    /// Default constructor.  
    AcdTkrHitPoca();
    
    /// Constructor for use in transient -> persistent conversion 
    /// Takes arguements as they are stored in ROOT
    AcdTkrHitPoca(const idents::AcdId& acdId, int trackIndex,
                  const float active2d[2], const float mips[2],
                  float vetoSigmaHit, float vetoSigmaProj, float vetoSigmaProp,
                  int volumePlane, float arcLengthToPlane, float cosTheta, 
                  const HepPoint3D& global, const float localPosition[2], 
                  const CLHEP::HepSymMatrix& localCovProj, const CLHEP::HepSymMatrix& localCovProp,
                  int volume, int region, float arcLength, 
                  float doca, float docaErrProj, float docaErrProp,
                  const Point& poca, const Vector& voca, const unsigned short flags[2]);

    /// Old Constructor for backwards compatiblity
    AcdTkrHitPoca( const idents::AcdId& acdId, int trackIndex, 
                   const Event::AcdTkrLocalCoords& local, const Event::AcdPocaData& pocaData );

    
    /// Copy constructor
    AcdTkrHitPoca(const Event::AcdTkrHitPoca& params);
    
    /// Assignment operator
    Event::AcdTkrHitPoca& operator=(const Event::AcdTkrHitPoca& params);

    /// Equality test operator, requires identity
    bool operator==(const Event::AcdTkrHitPoca& other) const {
      return this == &other;
    }

    /// Comparison operator, sorts by vetoSigma
    bool operator<(const Event::AcdTkrHitPoca& other) const;    

    /// Return the AcdId of the hit tile or ribbon
    inline const idents::AcdId& getId() const { return m_id; }
    
    /// Return the index of the associated track
    inline int trackIndex() const {
      return m_trackIndex;
    }

    /// Return the mips associated with PMT A
    inline float mipsPmtA() const { 
      return m_mips[0];
    }

    /// Return the mips associated with PMT B
    inline float mipsPmtB() const { 
      return m_mips[1];
    }

    /// Return the flags associated with PMT A
    inline unsigned short flagsPmtA() const { 
      return m_flags[0];
    }

    /// Return the flags associated with PMT B
    inline unsigned short flagsPmtB() const { 
      return m_flags[1];
    }

    /// Return true if more than 0.001 MIPs
    inline bool hasHit() const {
      return ( m_mips[0] > 0.001 || m_mips[1] > 0.001 );
    }

    /// combine the sigma from the hit with the sigma from the track
    float vetoSigma2() const;

    /// An estimator of the number of sigma needed for this track to be a true MIP signal
    inline float vetoSigmaHit() const {
      return m_vetoSigmaHit;
    }

    /// An estimator of the number of sigma needed for this track to hit this element
    inline float vetoSigmaProj() const {
      return m_vetoSigmaProj;
    }

    /// An estimator of the number of sigma needed for this track to hit this element
    inline float vetoSigmaProp() const {
      return m_vetoSigmaProp;
    }
        
    /// set all the values
    void set(const idents::AcdId& acdId, int trackIndex,
             const float mips[2],
             float vetoHit, float vetoProj, float vetoProp,
             const unsigned short flags[2]);                

    /// reset all the values to their default
    virtual void ini();
    
    /// Print out this structure
    virtual void writeOut(MsgStream& stream) const;
    
  private:
  
    /// The ID of the hit tile
    idents::AcdId m_id;
    
    /// The index of the associated track
    int m_trackIndex;
       
    /// The mip values associated with the two pmts
    float m_mips[2];

    /// The status bits from the ACD hit
    unsigned short m_flags[2];

    ///  An estimator of the number of sigma needed for this hit to be a true MIP signal
    float m_vetoSigmaHit;
    
    ///  An estimator of the number of sigma needed for this track to hit this element
    float m_vetoSigmaProj;

    ///  An estimator of the number of sigma needed for this track to hit this element
    float m_vetoSigmaProp;

  };


  /*! 
   * @class AcdTkrHitPocaCol
   *
   *  @brief TDS class a collection of AcdTkrHitPoca objects
   *  
   * It inherits from DataObject
   * (to be a part of Gaudi TDS) and from std::vector of pointers
   * to AcdTkrHitPoca objects. Some methods just rename the std::vector
   * methods with the same functionality for backward compartibility.
   *
   *
   * @author Eric Charles
   *
   * @todo replace this class by typedef to ObjectVector, make corresponding
   *       changes in AcdReconAlg
   */
    

  class AcdTkrHitPocaCol : public DataObject, public std::vector<AcdTkrHitPoca*> 
  {
  public:
    
    /// Default constructor.  Builds empty collection
    AcdTkrHitPocaCol() { clear();}
    
    /// "copy" constructor.  Take ownerships of a vector of AcdTkrHitPoca
    AcdTkrHitPocaCol(const std::vector<AcdTkrHitPoca*>& acdhits);
  
    /// destructor - deleting the hits pointed
    /// by the vector elements
    ~AcdTkrHitPocaCol() { delTkrHitPocas();}

    
    // GAUDI members to be use by the converters
    static const CLID& classID() {return CLID_AcdTkrHitPocaCol;}
    virtual const CLID& clID() const {return classID();}

    /// takes ownership of a vector AcdTkrHitPoca
    void init(std::vector<AcdTkrHitPoca*>& acdhits) {
      for ( std::vector<AcdTkrHitPoca*>::iterator itr = acdhits.begin(); itr != acdhits.end(); itr++ ) {
        push_back(*itr);
      }
    }
    
    /// add new AcdTkrHitPoca
    void add(AcdTkrHitPoca* cl) {push_back(cl);}
    
    /// get the number of pocas in collection
    int num()                  const {return size();}
    
    /// get pointer to the poca with a given number 
    AcdTkrHitPoca * getTkrHitPoca(int i) const {return operator[](i);}
    
    /// delete all pocas pointed by the vector elements
    void delTkrHitPocas();
    
    /// write information for all pocas to the ascii file 
    /// for debugging purposes
    virtual void writeOut(MsgStream& stream) const;
    
  protected:
    
    /// does the same function as clear() 
    virtual void ini();
        
  };


};
#endif
