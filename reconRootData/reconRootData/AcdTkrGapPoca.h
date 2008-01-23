#ifndef ROOT_ACDTKRGAPPOCA_H
#define ROOT_ACDTKRGAPPOCA_H

#include "TObject.h"
#include "AcdTkrLocalCoords.h"
#include "AcdPocaData.h"

/** 
 * @class AcdTkrGapPoca
 * @brief ROOT object which stores information about the Point of Closest Approach (POCA) between 
 * an extrapolated track and a Gap in the ACD
 *  
 * Most of the structure of the object comes from the base classes AcdTkrLocalCoords and AcdPocaData
 * 
 * The class adds only enough information to define the involved in the POCA.
 *    - const AcdGapId& getId()  
 *      - which returns the ID of the gap
 *    - const int getTrackIndex()  
 *      - which returns the index of the track which did the hitting
 *
 * 
 * @author Eric Charles
 *
  */

#include "commonRootData/idents/AcdGapId.h"

class AcdTkrGapPoca : public TObject, public AcdTkrLocalCoords, public AcdPocaData {

public:

  /// Default constructor.  
  AcdTkrGapPoca();
  
  /// Constructor for use in transient -> persistent conversion 
  /// Takes arguements as they are stored in ROOT
  AcdTkrGapPoca(const AcdGapId& gapId, int trackIndex,
		const AcdTkrLocalCoords& local,
		const AcdPocaData& pocaData);

  /// Copy constructor
  AcdTkrGapPoca(const AcdTkrGapPoca& params);

  /// Assignment operator
  AcdTkrGapPoca& operator=(const AcdTkrGapPoca& params);

  /// Return the AcdGapId of the Gap
  inline const AcdGapId& getGapId() const { return m_id; }
  
  /// Return the index of the associated track
  inline Int_t getTrackIndex() const {
    return m_trackIndex;
  }

  /// set all the values
  void initialize(const AcdGapId& acdId, int trackIndex,
		  const AcdTkrLocalCoords& local,
		  const AcdPocaData& pocaData);
  
  /// reset all the values to their default
  virtual void Clear(Option_t*);
  
  /// Print out this structure
  virtual void Print(Option_t *option) const;
  
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
  Bool_t CompareInRange( const AcdTkrGapPoca&, const std::string & name = "" ) const ; // for tests   

private:
  
  /// The ID of the Gap
  AcdGapId m_id;
  
  /// The index of the associated track
  Int_t m_trackIndex;

  ClassDef(AcdTkrGapPoca,1)
  
};
#endif
