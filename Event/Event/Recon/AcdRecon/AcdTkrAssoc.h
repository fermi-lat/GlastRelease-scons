#ifndef Event_ACDTKRASSOC_H
#define Event_ACDTKRASSOC_H

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

static const CLID& CLID_AcdTkrAssocCol = InterfaceID("AcdTkrAssocCol", 1, 0);

/**
 *  @class Event::AcdTkrAssoc
 *
 *  @brief
 *  \author Eric Charles
 *
 * $Header$
 **/

namespace Event
{

  class AcdTkrHitPoca;
  class AcdTkrGapPoca;
  class AcdTkrPoint;

  class AcdTkrAssoc : virtual public ContainedObject {
    
  public:

    /// Default constructor
    AcdTkrAssoc();

    /// Copy constructor
    AcdTkrAssoc(const AcdTkrAssoc& other);

    /// Constructor for use in reconstruction, 
    AcdTkrAssoc(int index, bool up, float energy, 
                const HepPoint3D& start, const HepVector3D& dir, float arcLength,
                const HepSymMatrix& covStart, const HepSymMatrix& covEnd,
                int tkrSSDVeto, float cornerDoca);
    
    /// Destructor is trivial
    virtual ~AcdTkrAssoc() {};
    
    /// Direct access to parameters
    inline int getTrackIndex() const { return m_index; }

    inline bool getUpward() const { return m_upward; }

    inline float getEnergy() const { return m_energy; }

    inline const HepPoint3D& getStart() const { return m_start; }

    inline const HepVector3D& getDir() const { return m_dir; }

    inline float getArcLength() const { return m_arcLength; }

    inline int getTkrSSDVeto() const { return m_tkrSSDVeto; }

    inline const HepSymMatrix& getCovStart() const { return m_cov_start; }

    inline const HepSymMatrix& getCovEnd() const { return m_cov_end; }

    inline unsigned nHitPoca() const { return m_hitPocae.size(); }

    inline unsigned nGapPoca() const { return m_gapPocae.size(); }

    inline float getCornerDoca() const { return m_cornerDoca; }

    inline const AcdTkrHitPoca* getHitPoca(unsigned i = 0) const { return m_hitPocae.size() > i ? m_hitPocae[i] : 0; }

    inline const AcdTkrGapPoca* getGapPoca(unsigned i = 0) const { return m_gapPocae.size() > i ? m_gapPocae[i] : 0; }

    inline const AcdTkrPoint* getPoint() const { return m_point; }

    /// set everything at once
    void set(int index, bool up, float energy, 
             const HepPoint3D& start, const HepVector3D& dir, float arcLength,
             const HepSymMatrix& covStart, const HepSymMatrix& covEnd,
             int tkrSSDVeto, float cornerDoca);

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

    HepSymMatrix    m_cov_start;

    HepSymMatrix    m_cov_end;     

    int             m_tkrSSDVeto;

    float           m_cornerDoca;

    std::vector<const AcdTkrHitPoca*> m_hitPocae;
    
    std::vector<const AcdTkrGapPoca*> m_gapPocae;

    AcdTkrPoint*    m_point;

  };

   
  /*! 
   * @class AcdTkrAssocCol
   *
   *  @brief TDS class to store the set of AcdTkrAssocs
   *  
   * It inherits from DataObject
   * (to be a part of Gaudi TDS) and from std::vector of pointers
   * to AcdTkrAssoc objects. Some methods just rename the std::vector
   * methods with the same functionality for backward compartibility.
   *
   * @author Eric Charles
   *
   * @todo replace this class by typedef to ObjectVector, make corresponding
   *       changes in AcdPha2MipTool
   */
    

  class AcdTkrAssocCol : public DataObject, public std::vector<AcdTkrAssoc*> 
  {
  public:

    /// Default constructor.  Builds empty collection
    AcdTkrAssocCol() { clear();}

    /// "copy" constructor.  Take ownerships of a vector of AcdTkrAssocs
    AcdTkrAssocCol(const std::vector<AcdTkrAssoc*>& acdhits);
    
    /// destructor - deleting the hits pointed
    /// by the vector elements
    ~AcdTkrAssocCol() { delTkrAssocs();}
            
    // GAUDI members to be use by the converters
    static const CLID& classID() {return CLID_AcdTkrAssocCol;}
    virtual const CLID& clID() const {return classID();}

    /// takes ownership of a vector AcdTkrAssoc
    void init(std::vector<AcdTkrAssoc*>& other) {
      for ( std::vector<AcdTkrAssoc*>::iterator itr = other.begin(); itr != other.end(); itr++ ) {
        push_back(*itr);
      }
    }
   
    /// Add a new hit
    void add(AcdTkrAssoc* cl) {push_back(cl);}
    
    /// get the number of hits in collection
    int num()                  const {return size();}
    
    /// get pointer to the hit with given number 
    AcdTkrAssoc * getTkrAssoc(int i) const {return operator[](i);}
    
    /// delete all Track associations pointed by the vector elements
    void delTkrAssocs();
    
    /// write information for all hits to the ascii file 
    /// for debugging purposes
    virtual void writeOut(MsgStream& stream) const;
    
  protected:
    
    /// does the same function as clear() 
    virtual void ini();
        
  };

}

#endif
