#ifndef Event_ACDTKRPOCA_H
#define Event_ACDTKRPOCA_H

#include <vector>

#include "Event/TopLevel/Definitions.h"

#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"

#include "idents/AcdId.h"
#include "geometry/Point.h"

#include "Event/Recon/TkrRecon/TkrTrackParams.h"

class MsgStream;

static const CLID& CLID_AcdTkrPocaCol = InterfaceID("AcdTkrPocaCol", 1, 0);

/**
*  @class Event::AcdTkrPoca
*
*  @brief Deprecated!!
*
*  \author Eric Charles
*
* $Header$
*/

namespace Event
{
  
  class AcdTkrPoca : virtual public ContainedObject {

  public:
    
    /// what region of a tile does the poca occur in
    typedef enum {
      NONE_TILE = 0,
      HIT_TILE = 1,
      MISS_TILE = 2,      
      HIT_TOP, HIT_RIGHT, HIT_BOTTOM, HIT_LEFT,
      MISS_TOP, MISS_TOP_RIGHT, MISS_RIGHT, MISS_BOTTOM_RIGHT,
      MISS_BOTTOM, MISS_BOTTOM_LEFT, MISS_LEFT, MISS_TOP_LEFT } DocaRegionTile;
    
    /// what region of a ribbon does the poca occur in
    typedef enum {
      NONE_RIBBON = 0,
      HIT_RIBBON = 1,
      MISS_RIBBON = 2,      
      MISS_POSITIVE, MISS_NEGATIVE } DocaRegionRibbon;

  public:
    
    /// Default constructor.  
    AcdTkrPoca();
    
    /// Constructor for use in persistent -> transient conversion and reconstruction
    /// Takes arguements as they are stored in ROOT and caluclated by AcdPocaTool
    AcdTkrPoca(const idents::AcdId& acdId, int trackIndex,
               double m_doca, double m_docaErr, unsigned docaRegion,
               const Point& poca, const Event::TkrTrackParams& paramsAtPoca);

    /// Destructor is trivial
    virtual ~AcdTkrPoca() {;}

    /// set all the values
    void set(const idents::AcdId& acdId, int trackIndex,
             double doca, double docaErr, unsigned docaRegion,
             const Point& poca, const Event::TkrTrackParams& paramsAtPoca);

    /// set only some of the values
    inline void setDocaErr(double val) { m_docaErr = val; };

    /// Return the AcdId of the hit tile or ribbon
    inline const idents::AcdId& acdId() const {
      return m_acdId;
    }

    /// Return the index of the associated track
    inline int trackIndex() const {
      return m_trackIdx;
    }

    /// Return the distance of closest approach (doca)
    /// This is calculated in 3D.  See AcdPocaTool for details
    inline float doca() const {
      return m_doca;
    }

    /// Return the error on distance of closest approach (doca)
    /// This is calulated the full Kalman Filter propagation of the track
    /// to the POCA and the projection of the propagated covarience matrix 
    /// along the line of the doca.
    inline float docaErr() const {
      return m_docaErr;
    }
    
    /// Return a code showing which region of the tile or ribbon the poca occured in. 
    inline unsigned docaRegion() const {
      return m_docaRegion;
    }

    /// Return the POCA (in global coordinates)
    const Point& poca() const {
      return m_poca;
    }

    /// Return the kalman propagated track parameters at the POCA
    const TkrTrackParams& paramsAtPoca() const {
      return m_paramsAtPoca;
    }

    /// Print out this structure on a stream
    virtual void writeOut(MsgStream& stream) const;
    
  protected:

    /// Reset all the values to their defaults
    virtual void ini();
    
  private:
    
    /// The tile or ribbon hit
    idents::AcdId m_acdId;

    /// The track index  
    int m_trackIdx;    

    /// The distance of closest approach (doca)
    float m_doca;

    /// The error on distance of closest approach (doca)
    float m_docaErr;
    
    /// a code showing which region of the tile or ribbon the poca occured in. 
    unsigned m_docaRegion;
    
    /// The point of closest approach of the track to the tile or ribbon
    Point m_poca;

    /// The track parameters (and covarience matrix) at the POCA
    Event::TkrTrackParams m_paramsAtPoca;
    
  };

   
  /*! 
   * @class AcdTkrPocaCol
   *
   *  @brief TDS class a collection of AcdTkrPoca objects
   *  
   * It inherits from DataObject
   * (to be a part of Gaudi TDS) and from std::vector of pointers
   * to AcdTkrPoca objects. Some methods just rename the std::vector
   * methods with the same functionality for backward compartibility.
   *
   *
   * @author Eric Charles
   *
   * @todo replace this class by typedef to ObjectVector, make corresponding
   *       changes in AcdReconAlg
   */
    

  class AcdTkrPocaCol : public DataObject, public std::vector<AcdTkrPoca*> 
  {
  public:
    
    /// Default constructor.  Builds empty collection
    AcdTkrPocaCol() { clear();}
    
    /// "copy" constructor.  Take ownerships of a vector of AcdTkrPoca
    AcdTkrPocaCol(const std::vector<AcdTkrPoca*>& acdhits);
  
    /// destructor - deleting the hits pointed
    /// by the vector elements
    ~AcdTkrPocaCol() { delTkrPocas();}
        
    
    // GAUDI members to be use by the converters
    static const CLID& classID() {return CLID_AcdTkrPocaCol;}
    virtual const CLID& clID() const {return classID();}
    
    /// add new AcdTkrPoca
    void add(AcdTkrPoca* cl) {push_back(cl);}
    
    /// get the number of pocas in collection
    int num()                  const {return size();}
    
    /// get pointer to the poca with a given number 
    AcdTkrPoca * getTkrPoca(int i) const {return operator[](i);}
    
    /// delete all pocas pointed by the vector elements
    void delTkrPocas();
    
    /// write information for all pocas to the ascii file 
    /// for debugging purposes
    virtual void writeOut(MsgStream& stream) const;
    
  protected:
    
    /// does the same function as clear() 
    virtual void ini();
        
  };

}

#endif
