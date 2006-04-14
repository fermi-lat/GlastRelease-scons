#ifndef ROOT_CCSDS_H
#define ROOT_CCSDS_H 1

#include "Riostream.h"
#include "TObject.h"

/** @class Ccsds 
* @brief Encapsulate CCSDS values 
*
*  The Ccsds consists of:
*  - scid
*  - apid.
*  - utc 
*
* $Header$
*/

class Ccsds : public TObject {
  
public:

  /// Default c'tor.  Assigns sentinel or null values to all fields
  Ccsds()
    :m_scid(0),m_apid(0),m_utc(0.0) {
  }
  
  /// Standard c'tor.  Takes input values for all fields
  Ccsds( int scid, int apid, double utc)
    :m_scid(scid),m_apid(apid),
     m_utc(utc) {
  }
 
  /// Copy c'tor.  Just copy all values.  
  Ccsds( const Ccsds& other )
    :TObject(other),
     m_scid(other.getScid()),
     m_apid(other.getApid()),
     m_utc(other.getUtc()) {
  }

  virtual ~Ccsds(){
  }
  
  inline Int_t getScid() const { return m_scid; }
  
  inline Int_t getApid() const { return m_apid; } 
  
  inline Double_t getUtc() const { return m_utc; }
  
  /// set everything at once
  inline void initialize(int scid, int apid, double utc) {
    m_scid = scid;
    m_apid = apid;
    m_utc = utc;
  }
  
  /// Reset function
  void Clear(Option_t* /* option="" */) {
    m_scid = 0;
    m_apid = 0;
    m_utc = 0.0;
  }

  /// ROOT print function
  void Print(Option_t* /* option="" */) const;
  
  /// Create a fake for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 
  
  /// Compare to another in tests
  Bool_t CompareInRange( const Ccsds&, const std::string & name = "" ) const ; 

private:
  
  Int_t m_scid;
  Int_t m_apid;
  Double_t m_utc;
 
  ClassDef(Ccsds,1) // CCSDS data

};

#endif    // ROOT_METAEVENT
