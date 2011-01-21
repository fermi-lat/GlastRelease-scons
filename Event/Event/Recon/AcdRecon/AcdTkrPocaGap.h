#ifndef Event_ACDTKRPOCAGAP_H
#define Event_ACDTKRPOCAGAP_H

#include <vector>

#include "Event/TopLevel/Definitions.h"

#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"

#include "idents/AcdId.h"
#include "geometry/Point.h"

#include "Event/Recon/TkrRecon/TkrTrackParams.h"

class MsgStream;

static const CLID& CLID_AcdTkrPocaGapCol = InterfaceID("AcdTkrPocaGapCol", 1, 0);

/**
*  @class Event::AcdTkrPocaGap
*
*  @brief Deprecated!!
*  
*  \author Eric Charles
*
* $Header$
*/

namespace Event
{
  
  class AcdTkrPocaGap : virtual public ContainedObject {

  public:
    
    /// Default constructor.  
    AcdTkrPocaGap();
    
    /// Constructor for use in persistent -> transient conversion and reconstruction
    /// Takes arguements as they are stored in ROOT and caluclated by AcdPocaGapTool
    AcdTkrPocaGap(const idents::AcdId& acdId, int trackIndex,
               double m_doca, double m_docaErr,
               const Point& PocaGap, const Vector& PocaGapVector, 
               const Event::AcdTkrLocalCoords& local);

    /// Destructor is trivial
    virtual ~AcdTkrPocaGap() {;}

    /// Return the AcdId of the hit tile or ribbon
    inline const idents::AcdId& acdId() const {
      return m_acdId;
    }

    /// Return the index of the associated track
    inline int trackIndex() const {
      return m_trackIdx;
    }

    /// Return the distance of closest approach (doca)
    /// This is calculated in 3D.  See AcdPocaGapTool for details
    inline float doca() const {
      return m_doca;
    }

    /// Return the error on distance of closest approach (doca)
    /// This is calulated the full Kalman Filter propagation of the track
    /// to the POCAGAP and the projection of the propagated covarience matrix 
    /// along the line of the doca.
    inline float docaErr() const {
      return m_docaErr;
    }
    
    /// Return the POCAGAP (in global coordinates)
    const Point& PocaGap() const {
      return m_PocaGap;
    }

    /// Return the vector from the POCAGAP to the edge of the detector element
    const Vector& PocaGapVector() const {
      return m_PocaGapVector;
    }

    /// set all the values at once
    void set(const idents::AcdId& acdId, int trackIndex,
             double m_doca, double m_docaErr,
             const Point& PocaGap, const Vector& PocaGapVector, 
             const Event::AcdTkrLocalCoords& local);
    
    // set the individual values (uncomment as needed)
 
    //inline void setAcdId(const idents::AcdId& val) { m_acdId = val; };
    //inline void setTrackIdx(int val) { m_trackIdx = val; };
    inline void setDocaErr(float val) { m_docaErr = val; };
    //inline void setPocaGap(const Point& val) { m_PocaGap = val; };
    //inline void setPocaGapVector(const Vector& val) { m_PocaGapVector = val; };
    

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
    
    /// The point of closest approach of the track to the gap
    Point m_PocaGap;

    /// The error on distance of closest approach (doca)
    float m_docaErr;

    /// The vector from the POCA to the gap
    Vector m_PocaGapVector;

    /// The 2D Local coords (in the plane of the elements)
    Event::AcdTkrLocalCoords m_local;
    
  };

   
  /*! 
   * @class AcdTkrPocaGapCol
   *
   *  @brief TDS class  to store the results of the reconstruction performed
   *  
   * It inherits from DataObject
   * (to be a part of Gaudi TDS) and from std::vector of pointers
   * to AcdTkrPocaGap objects. Some methods just rename the std::vector
   * methods with the same functionality for backward compartibility.
   *
   *
   * @author Eric Charles
   *
   * @todo replace this class by typedef to ObjectVector, make corresponding
   *       changes in AcdTrkIntersectionAlg
   */
    

  class AcdTkrPocaGapCol : public DataObject, public std::vector<AcdTkrPocaGap*> 
  {
  public:
    
    /// Default constructor.  Builds empty collection
    AcdTkrPocaGapCol() { clear();}
    
    /// "copy" constructor.  Take ownerships of a vector of AcdTkrPocaGap
    AcdTkrPocaGapCol(const std::vector<AcdTkrPocaGap*>& acdhits);
  
    /// destructor - deleting the hits pointed
    /// by the vector elements
    ~AcdTkrPocaGapCol() { delTkrPocaGaps();}
        
    
    // GAUDI members to be use by the converters
    static const CLID& classID() {return CLID_AcdTkrPocaGapCol;}
    virtual const CLID& clID() const {return classID();}
    
    /// add new AcdTkrPocaGap
    void add(AcdTkrPocaGap* cl) {push_back(cl);}
    
    /// get the number of PocaGaps in collection
    int num()                  const {return size();}
    
    /// get pointer to the PocaGap with a given number 
    AcdTkrPocaGap * getTkrPocaGap(int i) const {return operator[](i);}
    
    /// delete all PocaGaps pointed by the vector elements
    void delTkrPocaGaps();
    
    /// write information for all PocaGaps to the ascii file 
    /// for debugging purposes
    virtual void writeOut(MsgStream& stream) const;
    
  protected:
    
    /// does the same function as clear() 
    virtual void ini();
        
  };

}

#endif
