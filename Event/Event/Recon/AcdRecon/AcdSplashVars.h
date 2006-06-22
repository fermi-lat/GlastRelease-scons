#ifndef EVENT_ACDBEAMVARS_H
#define EVENT_ACDBEAMVARS_H

#include <vector>
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"

#include "idents/AcdId.h"

#include "geometry/Point.h"
#include "geometry/Vector.h"

class MsgStream;

/** 
 * @class AcdSplashVars
 * @brief Root object information about the Point of Closest Approach (POCA) between an extrapolated track
 *  and a hit Acd element (tile or ribbon).  This POCA is calculated in 3D.  The doca is defined to be positive 
 *  if the track goes inside the active distance and negative otherwise
 *  
 * This class should be a duplicate of Event::AcdSplashVars
 * 
 * @author Eric Charles
 *
 * $Header$
 */

static const CLID& CLID_AcdSplashVarsCol = InterfaceID("AcdSplashVarsCol", 1, 0);

namespace Event {

  class AcdSplashVars {

  public:
    
    /// Default constructor.  
    AcdSplashVars();
    
    /// Constructor for use in transient -> persistent conversion 
    /// Takes arguements as they are stored in ROOT
    AcdSplashVars(const idents::AcdId& acdId, int trackIndex, 
		const Point& calEntryPoint, const Vector& calEntryVector,
		const float& tileSolidAngle, const float& weightedTrackAngle,
		const float& weightedPathlength);
    
    /// Copy constructor
    AcdSplashVars(const AcdSplashVars& params);
    
    /// Assignment operator
    AcdSplashVars& operator=(const AcdSplashVars& params);

    /// D'tor
    virtual ~AcdSplashVars() {;}
    
    /// Return the AcdId of the hit tile or ribbon
    inline const idents::AcdId& getId() const { return m_id; }
    
    /// Return the index of the associated track
    inline int getTrackIndex() const {
      return m_trackIndex;
    }

    /// return the Vector from the point the track enters the calorimeter to the tile center
    inline const Point& calEntryPoint() const { return m_calEntryPoint; }
    
    /// return the track vector at the point the track enters the calorimeter
    inline const Vector& calEntryVector() const { return m_calEntryVector; }
    
    /// return total solid angle of the tile, seen from the track entry point
    inline const float&  tileSolidAngle() const { return m_tileSolidAngle; }
    
    /// return the average of the angle between the reconstructed track and the line         
    /// connecting the track entry point in the calorimeter at that point
    /// (weighted by the solid angle of the element)                         
    inline const float&  weightedTrackAngle() const { return m_weightedTrackAngle; }
    
    /// return the average of the pathlength inside the tile along the path from the 
    /// track entry point in the calorimeter to the tile
    /// (weighted by the solid angle of the element)     
    inline const float&  weightedPathlength() const { return m_weightedPathlength; }
    
    /// set all the values
    void set(const idents::AcdId& acdId, int trackIndex, 
	     const Point& calEntryPoint, const Vector& calEntryVector,
	     const float& tileSolidAngle, const float& weightedTrackAngle,
	     const float& weightedPathlength);
    
    /// reset all the values to their default
    virtual void ini();
    
    /// Print out this structure
    virtual void writeOut(MsgStream& stream) const;
    
  private:
    
    /// The ID of the hit tile
    idents::AcdId m_id;
    
    /// The index of the associated track
    int m_trackIndex;
    
    /// The Vector from the point the track enters the calorimeter to the tile center
    Point m_calEntryPoint;
    
    /// The track vector at the point the track enters the calorimeter
    Vector m_calEntryVector;
    
    /// The total solid angle of the tile, seen from the track entry point
    float m_tileSolidAngle;
    
    /// The average of the angle between the reconstructed track and the line         
    /// connecting the track entry point in the calorimeter at that point
    /// (weighted by the solid angle of the element)                         
    float m_weightedTrackAngle;
    
    /// the average of the pathlength inside the tile along the path from the 
    /// track entry point in the calorimeter to the tile
    /// (weighted by the solid angle of the element)     
    float m_weightedPathlength;
  
  };

  class AcdSplashVarsCol : public DataObject, public std::vector<AcdSplashVars*> 
  {
  public:
    
    /// Default constructor.  Builds empty collection
    AcdSplashVarsCol() { clear();}
    
    /// "copy" constructor.  Take ownerships of a vector of AcdSplashVars
    AcdSplashVarsCol(const std::vector<AcdSplashVars*>& vars);
  
    /// destructor - deleting the hits pointed
    /// by the vector elements
    ~AcdSplashVarsCol() { del();}
        
    
    // GAUDI members to be use by the converters
    static const CLID& classID() {return CLID_AcdSplashVarsCol;}
    virtual const CLID& clID() const {return classID();}
    
    /// add new AcdSplashVars
    void add(AcdSplashVars* cl) {push_back(cl);}
    
    /// get the number of pocas in collection
    int num()                  const {return size();}
    
    /// get pointer to the poca with a given number 
    AcdSplashVars * get(int i) const {return operator[](i);}
    
    /// delete all pocas pointed by the vector elements
    void del();
    
    /// write information for all pocas to the ascii file 
    /// for debugging purposes
    virtual void writeOut(MsgStream& stream) const;
    
  protected:
    
    /// does the same function as clear() 
    virtual void ini();
        
  };

    
};

#endif
