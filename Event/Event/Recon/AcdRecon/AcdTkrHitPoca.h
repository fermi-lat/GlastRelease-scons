#ifndef EVENT_ACDTKRHITPOCA_H
#define EVENT_ACDTKRHITPOCA_H

#include "AcdTkrLocalData.h"
#include "AcdPocaData.h"

/** 
 * @class AcdTkrHitPoca
 * @brief TDS object information about the Point of Closest Approach (POCA) between an extrapolated track
 *  and a hit Acd element (tile or ribbon).  This POCA is calculated in 3D.  The doca is defined to be positive 
 *  if the track goes inside the active distance and negative otherwise
 *  
 * This class should be a duplicate of reconRootData/AcdTkrHitPoca
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
		  const Event::AcdTkrLocalCoords& local,
		  const Event::AcdPocaData& pocaData);
    
    /// Copy constructor
    AcdTkrHitPoca(const Event::AcdTkrHitPoca& params);
    
    /// Assignment operator
    Event::AcdTkrHitPoca& operator=(const Event::AcdTkrHitPoca& params);

    /// Return the AcdId of the hit tile or ribbon
    inline const idents::AcdId& getId() const { return m_id; }
    
    /// Return the index of the associated track
    inline int trackIndex() const {
      return m_trackIndex;
    }
    
    /// set all the values
    void set(const idents::AcdId& acdId, int trackIndex,
	     const Event::AcdTkrLocalCoords& local,
	     const Event::AcdPocaData& pocaData);
    
    /// reset all the values to their default
    virtual void ini();
    
    /// Print out this structure
    virtual void writeOut(MsgStream& stream) const;
    
  private:
  
    /// The ID of the hit tile
    idents::AcdId m_id;
    
    /// The index of the associated track
    int m_trackIndex;
    
  };


  /*! 
   * @class AcdTkrHitPocaCol
   *
   *  @brief TDS class  to store the results of the reconstruction performed
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
   *       changes in AcdTrkIntersectionAlg
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
