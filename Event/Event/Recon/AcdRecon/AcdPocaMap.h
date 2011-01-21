#ifndef Event_ACDPOCAMAP_H
#define Event_ACDPOCAMAP_H

#include <set>
#include <map>

#include "Event/TopLevel/Definitions.h"
#include "idents/AcdId.h"

#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"

#include "Event/Recon/AcdRecon/AcdTkrHitPoca.h"

class MsgStream;

static const CLID& CLID_AcdPocaMap = InterfaceID("AcdPocaMap", 1, 0);

/**
*  @class Event::AcdPocaMap
*
*  @brief This class stores the relations between TkrTracks, AcdIds and AcdTkrHitPocas.
*
*  Associations are stored two ways:
*  - A map from track index to sets of AcdTkrHitPoca sorted by doca value (largest active distance first)
*  - A map from AcdId to sets of AcdTkrHitPoca sorted by doca value (largest active distance first)
*
*  The main various functions are:
*  - const Event::AcdPocaSet& getPocas( const idents::AcdId& acdId )
*    - get all the pocas associated with a single tile
*  - const Event::AcdTkrHitPoca* getBestPoca(const idents::AcdId& acdId) 
*    - get the best poca associated with a single tile (can return null)
*  - const Event::AcdPocaSet& getPocas(int track)
*    - get all the pocas associated with a single track
*  - const Event::AcdTkrHitPoca* getBestPoca(int track) 
*    - get the best poca associated with a single track (can return null)
*  - const Event::AcdTkrHitPoca* getPoca(const idents::AcdId& acdId, int track)
*    - get the poca between a track and a tile (can return null)
*
*  
*  \author Eric Charles
*
* $Header$ */

namespace Event
{
  
  /// This allows us to sort sets of AcdTkrHitPoca by the doca value
  struct AcdPocaMore {
    bool operator()(const AcdTkrHitPoca* ptr1, const AcdTkrHitPoca* ptr2) const {        
      return ptr1->getDoca() > ptr2->getDoca();
    }
  };   
  
  //
  // Define an AcdPocaSet as an stl set, sorted by doca
  // The sorting is done so the large value occurs first.
  typedef std::set<const Event::AcdTkrHitPoca*, Event::AcdPocaMore> AcdPocaSet;


  
  class AcdPocaMap : virtual public DataObject {

  public:
    
    /// default constructor, builds any empty map
    AcdPocaMap();

    /// clear out the map
    void clear() {
      ini();
    }

    /// adds a single poca to the map
    void add(const Event::AcdTkrHitPoca& poca);

    /// get all the pocas associated with a single tile
    const Event::AcdPocaSet& getPocas(const idents::AcdId& acdId) const;
    /// get the best poca associated with a single tile (can return null)
    const Event::AcdTkrHitPoca* getBestPoca(const idents::AcdId& acdId) const;

    /// get all the pocas associated with a single track
    const Event::AcdPocaSet& getPocas(int track) const;
    /// get the best poca associated with a single track (can return null)
    const Event::AcdTkrHitPoca* getBestPoca(int track) const;

    /// get the poca between a track and a tile (can return null)
    const Event::AcdTkrHitPoca* getPoca(const idents::AcdId& acdId, int track) const;

  protected:

    /// clear the map
    virtual void ini();

    /// Print out this structure on a stream
    virtual void writeOut(MsgStream& stream) const;    
    
  private:
    
    /// map from tiles to pocas
    std::map< idents::AcdId, Event::AcdPocaSet > m_tileToPocaMap;
    
    /// map from tracks to pocas
    std::map< int, Event::AcdPocaSet > m_trackToPocaMap;
    
  };

}

#endif
