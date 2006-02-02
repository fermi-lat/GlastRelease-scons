#ifndef ROOT_DATAGRAMINFO_H
#define ROOT_DATAGRAMINFO_H 1

#include <iostream>
#include "TObject.h"

#include "enums/Lsf.h"

/** @class DatagramInfo
* @brief encapsulate the Datagram id parts of the event context
*
* 
*
* $Header$
*/

    
class DatagramInfo : public TObject {
    
public:
  
  DatagramInfo( )
    : m_openAction(enums::Lsf::Open::Unspecified), m_openReason(enums::Lsf::Open::Unknown), 
      m_crate(enums::Lsf::NoCrate), m_mode(enums::Lsf::NoMode), 
      m_closeAction(enums::Lsf::Close::Unspecified), m_closeReason(enums::Lsf::Close::Unknown), 
      m_datagrams(0),m_modeChanges(0){
  };
  
  DatagramInfo( enums::Lsf::Open::Action oAction, enums::Lsf::Open::Reason oReason, 
		enums::Lsf::Crate c, enums::Lsf::Mode m, 
		enums::Lsf::Close::Action cAction, enums::Lsf::Close::Reason cReason, 
		UInt_t datagrams, UInt_t modeChanges)
    : m_openAction(oAction), m_openReason(oReason), 
      m_crate(c), m_mode(m), 
      m_closeAction(cAction), m_closeReason(cReason),
      m_datagrams(datagrams),m_modeChanges(modeChanges){
  }
  
  ~DatagramInfo() {
  }
  
  DatagramInfo( const DatagramInfo& other ) 
    : TObject(other),
    m_openAction(other.openAction()), m_openReason(other.openReason()),
    m_crate(other.crate()), m_mode(other.mode()), 
    m_closeAction(other.closeAction()), m_closeReason(other.closeReason()),
    m_datagrams(other.datagrams()),m_modeChanges(other.modeChanges()){
  }
  
  inline DatagramInfo& operator=(const DatagramInfo& other) {
    initialize(other.openAction(),other.openReason(),
	       other.crate(),other.mode(),
	       other.closeAction(),other.closeReason(),
	       other.datagrams(),other.modeChanges());
    return *this;
  }
  
  /// number of mode changes since the start of the run
  inline UInt_t modeChanges() const { return m_modeChanges; };
  
  /// number of datagrams sent since the start of the run
  /// this is identical to the datagram sequence number
  inline UInt_t datagrams() const { return m_datagrams; };
  
  /// The action that caused the datagram to be opened
  inline enums::Lsf::Open::Action openAction() const { return m_openAction; }
  
  ///The reason this datagram was opened
  inline enums::Lsf::Open::Reason openReason() const { return m_openReason; }
  
  /// Source that this datagram came from
  inline enums::Lsf::Crate crate() const { return m_crate; }
  
  /// Operating mode the LAT was in when the data for this data were acquired
  inline enums::Lsf::Mode mode() const { return m_mode; } 
  
  /// The action that caused the datagram to be closed
  inline enums::Lsf::Close::Action closeAction() const { return m_closeAction; }
  
  ///The reason this datagram was closed
  inline enums::Lsf::Close::Reason closeReason() const { return m_closeReason; }
  
  /// set everything at once
  inline void initialize(enums::Lsf::Open::Action oAction, enums::Lsf::Open::Reason oReason,
			 enums::Lsf::Crate c, enums::Lsf::Mode m, 
			 enums::Lsf::Close::Action cAction, enums::Lsf::Close::Reason cReason,
			 UInt_t datagrams, UInt_t modeChanges) {
    m_openAction = oAction;
    m_openReason = oReason;
    m_crate = c;
    m_mode = m;
    m_closeAction = cAction;
    m_closeReason = cReason;
    m_datagrams = datagrams;
    m_modeChanges = modeChanges;
  }
  
  // set the individual data members
  inline void setModeChanges( UInt_t val ) { m_modeChanges = val; };
  inline void setDatagrams( UInt_t val ) { m_datagrams = val; }; 
  inline void setOpenAction( enums::Lsf::Open::Action val ) { m_openAction = val; };
  inline void setOpenReason( enums::Lsf::Open::Reason val ) { m_openReason = val; };
  inline void setCrate( enums::Lsf::Crate val ) { m_crate = val; };
  inline void setMode( enums::Lsf::Mode val ) { m_mode = val; }; 
  inline void setCloseAction( enums::Lsf::Close::Action val ) { m_closeAction = val; };
  inline void setCloseReason( enums::Lsf::Close::Reason val ) { m_closeReason = val; };

  /// Reset function
  void Clear(Option_t* /* option="" */) {
     m_openAction = enums::Lsf::Open::Unspecified; 
     m_openReason = enums::Lsf::Open::Unknown;
     m_crate = enums::Lsf::NoCrate; 
     m_mode = enums::Lsf::NoMode;
     m_closeAction = enums::Lsf::Close::Unspecified; 
     m_closeReason = enums::Lsf::Close::Unknown; 
     m_datagrams = 0;
     m_modeChanges = 0;
  }

  /// ROOT print function
  void Print(Option_t* /* option="" */) const;

  /// Create a fake for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 

  /// Compare to another in tests
  Bool_t CompareInRange( const DatagramInfo&, const std::string & name = "" ) const ; 

private:
  
  enums::Lsf::Open::Action m_openAction;
  enums::Lsf::Open::Reason m_openReason;
  enums::Lsf::Crate m_crate;
  enums::Lsf::Mode m_mode;    
  
  enums::Lsf::Close::Action m_closeAction;
  enums::Lsf::Close::Reason m_closeReason;
  
  UInt_t m_datagrams;    
  UInt_t m_modeChanges; 

  ClassDef(DatagramInfo,1);

};

#endif    // ROOT_DATAGRAMINFO_H
