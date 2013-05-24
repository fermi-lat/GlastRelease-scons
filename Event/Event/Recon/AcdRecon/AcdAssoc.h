#ifndef Event_ACDASSOC_H
#define Event_ACDASSOC_H

#include <vector>

#include "Event/TopLevel/Definitions.h"

#include "idents/AcdId.h"


#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"

#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"


class MsgStream;

static const CLID& CLID_AcdAssocCol = InterfaceID("AcdAssocCol", 1, 0);
static const CLID& CLID_AcdTkrAssocCol = InterfaceID("AcdTkrAssocCol", 1, 0);
static const CLID& CLID_AcdCalAssocCol = InterfaceID("AcdCalAssocCol", 1, 0);

typedef HepGeom::Point3D<double> HepPoint3D;
typedef HepGeom::Vector3D<double> HepVector3D;

/**
 *  @class Event::AcdAssoc
 *
 *  @brief
 *  \author Eric Charles
 *  \author Alex Drlica-Wagner
 *
 * $Header$
 **/

namespace Event
{

  class AcdTkrHitPoca;
  class AcdTkrGapPoca;
  class AcdTkrPoint;

  class AcdAssoc : virtual public ContainedObject {
    
  public:

    /// Default constructor
    AcdAssoc();

    /// Copy constructor
    AcdAssoc(const AcdAssoc& other);

    /// Constructor for use in reconstruction, 
    AcdAssoc(int index, bool up, float energy, 
	     const HepPoint3D& start, const HepVector3D& dir, float arcLength,
	     const CLHEP::HepSymMatrix& covStart, const CLHEP::HepSymMatrix& covEnd,
	     int tkrSSDVeto, float cornerDoca,
	     float energy15, float energy30, float energy45,
	     float triggerEnergy15, float triggerEnergy30, float triggerEnergy45);
    
    /// Destructor is trivial
    virtual ~AcdAssoc() {};
    
    /// Direct access to parameters
    inline int getTrackIndex() const { return m_index; }

    inline bool getUpward() const { return m_upward; }

    inline float getEnergy() const { return m_energy; }

    inline const HepPoint3D& getStart() const { return m_start; }

    inline const HepVector3D& getDir() const { return m_dir; }

    inline float getArcLength() const { return m_arcLength; }

    inline int getTkrSSDVeto() const { return m_tkrSSDVeto; }

    inline const CLHEP::HepSymMatrix& getCovStart() const { return m_cov_start; }

    inline const CLHEP::HepSymMatrix& getCovEnd() const { return m_cov_end; }

    inline unsigned nHitPoca() const { return m_hitPocae.size(); }

    inline unsigned nGapPoca() const { return m_gapPocae.size(); }

    inline float getCornerDoca() const { return m_cornerDoca; }

    inline const AcdTkrHitPoca* getHitPoca(unsigned i = 0) const { return m_hitPocae.size() > i ? m_hitPocae[i] : 0; }

    inline const AcdTkrGapPoca* getGapPoca(unsigned i = 0) const { return m_gapPocae.size() > i ? m_gapPocae[i] : 0; }

    inline const AcdTkrPoint* getPoint() const { return m_point; }

    inline float GetEnergy15() const { return m_energy15; }
    
    inline float GetEnergy30() const { return m_energy30; }
    
    inline float GetEnergy45() const { return m_energy45; }
    
    inline float GetTriggerEnergy15() const { return m_triggerEnergy15; }
    
    inline float GetTriggerEnergy30() const { return m_triggerEnergy30; }
    
    inline float GetTriggerEnergy45() const { return m_triggerEnergy45; }  
    

    /// set everything at once
    void set(int index, bool up, float energy, 
             const HepPoint3D& start, const HepVector3D& dir, float arcLength,
             const CLHEP::HepSymMatrix& covStart, const CLHEP::HepSymMatrix& covEnd,
             int tkrSSDVeto, float cornerDoca,
	     float energy15, float energy30, float energy45,
	     float triggerEnergy15, float triggerEnergy30, float triggerEnergy45);

    /// add a hitPoca
    inline void addHitPoca(const AcdTkrHitPoca& poca) {
      m_hitPocae.push_back(&poca);
    }

    /// add a gapPoca
    inline void addGapPoca(const AcdTkrGapPoca& poca) {
      m_gapPocae.push_back(&poca);
    }    

    /// add a gapPoca
    inline void setPoint(AcdTkrPoint& point) {
      m_point = &point;
    }    

    /// Print out this structure on a stream
    virtual void writeOut(MsgStream& stream) const;
    
  protected:

    /// Reset all the values to their defaults
    virtual void ini();
    
  private:    

    int             m_index;
    
    bool            m_upward;

    float           m_energy;
    
    HepPoint3D      m_start;

    HepVector3D     m_dir;
 
    float           m_arcLength;

    CLHEP::HepSymMatrix    m_cov_start;

    CLHEP::HepSymMatrix    m_cov_end;     

    int             m_tkrSSDVeto;

    float           m_cornerDoca;

    std::vector<const AcdTkrHitPoca*> m_hitPocae;
    
    std::vector<const AcdTkrGapPoca*> m_gapPocae;

    AcdTkrPoint*    m_point;

    float           m_energy15;

    float           m_energy30;
    
    float           m_energy45;
    
    float           m_triggerEnergy15;
    
    float           m_triggerEnergy30;
    
    float           m_triggerEnergy45;

  };

   
  /*! 
   * @class AcdAssocCol
   *
   *  @brief TDS class to store the set of AcdAssocs
   *  
   * It inherits from DataObject
   * (to be a part of Gaudi TDS) and from std::vector of pointers
   * to AcdAssoc objects. Some methods just rename the std::vector
   * methods with the same functionality for backward compartibility.
   *
   * @author Eric Charles
   *
   * @todo replace this class by typedef to ObjectVector, make corresponding
   *       changes in AcdPha2MipTool
   */
    

  class AcdAssocCol : public DataObject, public std::vector<AcdAssoc*> 
  {
  public:

    /// Default constructor.  Builds empty collection
    AcdAssocCol() { clear();}

    /// "copy" constructor.  Take ownerships of a vector of AcdAssocs
    AcdAssocCol(const std::vector<AcdAssoc*>& acdhits);
    
    /// destructor - deleting the hits pointed
    /// by the vector elements
    ~AcdAssocCol() { delAssocs();}
            
    // GAUDI members to be use by the converters
    static const CLID& classID() {return CLID_AcdAssocCol;}
    virtual const CLID& clID() const {return classID();}

    /// takes ownership of a vector AcdAssoc
    void init(std::vector<AcdAssoc*>& other) {
      for ( std::vector<AcdAssoc*>::iterator itr = other.begin(); itr != other.end(); itr++ ) {
        push_back(*itr);
      }
    }
   
    /// Add a new association
    void add(AcdAssoc* cl) {push_back(cl);}
    
    /// get the number of hits in collection
    int num()                  const {return size();}
    
    /// get pointer to the hit with given number 
    AcdAssoc * getAssoc(int i) const {return operator[](i);}
    
    /// delete all associations pointed by the vector elements
    void delAssocs();
    
    /// write information for all hits to the ascii file 
    /// for debugging purposes
    virtual void writeOut(MsgStream& stream) const;
    
  protected:
    
    /// does the same function as clear() 
    virtual void ini();
        
  };

  class AcdTkrAssocCol : public AcdAssocCol
  {
  public:

    /// Default constructor.  Builds empty collection
    AcdTkrAssocCol() { clear();}

    /// "copy" constructor.  Take ownerships of a vector of AcdAssocs
    AcdTkrAssocCol(const std::vector<AcdAssoc*>& acdhits) : AcdAssocCol( acdhits ) {;}
    
    /// destructor - deleting the hits pointed
    /// by the vector elements
    ~AcdTkrAssocCol() { delTkrAssocs();}

    // GAUDI members to be use by the converters
    static const CLID& classID() {return CLID_AcdTkrAssocCol;}
    virtual const CLID& clID() const {return classID();}

    /// get pointer to the hit with given number 
    AcdAssoc * getTkrAssoc(int i) const {return operator[](i);}
    
    /// delete all Track associations pointed by the vector elements
    void delTkrAssocs() { delAssocs(); }
  };
  
  class AcdCalAssocCol : public AcdAssocCol
  {
  public:

    /// Default constructor.  Builds empty collection
    AcdCalAssocCol() { clear();}

    /// "copy" constructor.  Take ownerships of a vector of AcdAssocs
    AcdCalAssocCol(const std::vector<AcdAssoc*>& acdhits) : AcdAssocCol( acdhits ) {;}
    
    /// destructor - deleting the hits pointed
    /// by the vector elements
    ~AcdCalAssocCol() { delCalAssocs();}

    // GAUDI members to be use by the converters
    static const CLID& classID() {return CLID_AcdCalAssocCol;}
    virtual const CLID& clID() const {return classID();}

    /// get pointer to the hit with given number 
    AcdAssoc * getCalAssoc(int i) const {return operator[](i);}
    
    /// delete all Cal associations pointed by the vector elements
    void delCalAssocs() { delAssocs(); }
  };

}
#endif
