#ifndef ROOT_ACDRECONV2_H
#define ROOT_ACDRECONV2_H

#include "TObject.h"
#include "commonRootData/idents/AcdId.h"

#include "TClonesArray.h"

#include "AcdEventTopology.h"
#include "AcdHit.h"
#include "AcdTkrAssoc.h"

#include <iostream>


/** 
* @class AcdRecon
* @brief Root container for the reconstructed ACD information 
* 
* This class contains:

* @author Eric Charles
*
* $Header$
*/

class AcdReconV2 : public TObject {

public:

  AcdReconV2();
  
  virtual ~AcdReconV2();
  
  // Access to the data
  inline AcdEventTopology& getEventTopology() { return m_topology; }
  inline const AcdEventTopology& getEventTopology() const { return m_topology; }
  
  inline TClonesArray& getHitCol() { return *m_acdHitCol; }
  inline const TClonesArray& getHitCol() const { return *m_acdHitCol; }
  
  inline TClonesArray& getTkrAssocCol() { return *m_acdTkrAssocCol; }
  inline const TClonesArray& getTkrAssocCol() const { return *m_acdTkrAssocCol; }
  
  AcdHit* addAcdHit(const AcdHit& toAdd);
  AcdTkrAssoc* addAcdTkrAssoc(const AcdTkrAssoc& toAdd);

  void Clear(Option_t *option="");
  void Fake( Int_t ievent, Float_t randNum ) ; // for tests
  Bool_t CompareInRange( const AcdReconV2 &, const std::string & name = "" ) const ; // for tests
  void Print(Option_t *option="") const;
    
private:

  // Store the event topology
  AcdEventTopology m_topology;

  // Store the calibrated hits
  TClonesArray *m_acdHitCol;    
  
  // Store the track-tile( or ribbon) assocications
  TClonesArray *m_acdTkrAssocCol;

  ClassDef(AcdReconV2,1) // Acd Reconstruction data
};

#endif
