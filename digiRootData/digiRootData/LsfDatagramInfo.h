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
  

public:
  
  DatagramInfo( )
    : m_openReason(enums::Lsf::NoOpenReason), m_openRequester(enums::Lsf::NoOpenRequest), 
      m_crate(enums::Lsf::NoCrate), m_mode(enums::Lsf::NoMode), 
      m_closeReason(enums::Lsf::NoCloseReason), m_closeRequester(enums::Lsf::NoCloseRequest),
      m_datagrams(0),m_modeChanges(0){
  }
  
  DatagramInfo( enums::Lsf::OpenReason oReason, enums::Lsf::OpenRequester oReq, 
		enums::Lsf::Crate c, enums::Lsf::Mode m, 
		enums::Lsf::CloseReason cReason, enums::Lsf::CloseRequester cReq,
		UInt_t datagrams, UInt_t modeChanges)
    : m_openReason(oReason), m_openRequester(oReq), 
      m_crate(c), m_mode(m), 
      m_closeReason(cReason), m_closeRequester(cReq),
      m_datagrams(datagrams),m_modeChanges(modeChanges){
  }
  
  ~DatagramInfo() {
  }
  
  DatagramInfo( const DatagramInfo& other ) 
    : TObject(other),
      m_openReason(other.openReason()), m_openRequester(other.openRequester()),
      m_crate(other.crate()), m_mode(other.mode()), 
      m_closeReason(other.closeReason()), m_closeRequester(other.closeRequester()),
      m_datagrams(other.datagrams()),m_modeChanges(other.modeChanges()){
  }
  
  inline DatagramInfo& operator=(const DatagramInfo& other) {
    initialize(other.openReason(),other.openRequester(),
	       other.crate(),other.mode(),
	       other.closeReason(),other.closeRequester(),
	       other.datagrams(),other.modeChanges());
    return *this;
  }
  
  /// number of mode changes since the start of the run
  inline UInt_t modeChanges() const { return m_modeChanges; };
  
  /// number of datagrams sent since the start of the run
  /// this is identical to the datagram sequence number
  inline UInt_t datagrams() const { return m_datagrams; };
  
  /// Reason that this datagram was opened
  /// FIXME: should this be LsfEvent::DatagramInfo::Start for the first datagram in a run? 
  inline enums::Lsf::OpenReason openReason() const { return m_openReason; }
  
  /// Method used to open this datagram
  inline enums::Lsf::OpenRequester openRequester() const { return m_openRequester; }
  
  /// Source that this datagram came from
  inline enums::Lsf::Crate crate() const { return m_crate; }
  
    /// Operating mode the LAT was in when the data for this data were acquired
  inline enums::Lsf::Mode mode() const { return m_mode; } 
  
  /// Reason that this datagram was opened
  /// FIXME: should this be LsfEvent::DatagramInfo::Stop for the last datagram in a run?
  inline enums::Lsf::CloseReason closeReason() const { return m_closeReason; }
  
  /// Method used to close this datagram
  inline enums::Lsf::CloseRequester closeRequester() const { return m_closeRequester; }
  
  /// set everything at once
  inline void initialize(enums::Lsf::OpenReason oReason, enums::Lsf::OpenRequester oReq, 
			 enums::Lsf::Crate c, enums::Lsf::Mode m, 
			 enums::Lsf::CloseReason cReason, enums::Lsf::CloseRequester cReq,
			 UInt_t datagrams, UInt_t modeChanges) {
    m_openReason = oReason;
    m_openRequester = oReq;
    m_crate = c;
    m_mode = m;
    m_closeReason = cReason;
    m_closeRequester = cReq;
    m_datagrams = datagrams;
    m_modeChanges = modeChanges;
  }
  
  // set the individual data members
  inline void setModeChanges( UInt_t val ) { m_modeChanges = val; };
  inline void setDatagrams( UInt_t val ) { m_datagrams = val; }; 
  inline void setOpenReason( enums::Lsf::OpenReason val ) { m_openReason = val; };
  inline void setOpenRequester( enums::Lsf::OpenRequester val ) { m_openRequester = val; };
  inline void setCrate( enums::Lsf::Crate val ) { m_crate = val; };
  inline void setMode( enums::Lsf::Mode val ) { m_mode = val; }; 
  inline void setCloseReason( enums::Lsf::CloseReason val ) { m_closeReason = val; };
  inline void setCloseRequester( enums::Lsf::CloseRequester val ) { m_closeRequester = val; };
  
  /// Reset function
  void Clear(Option_t* /* option="" */) {
    m_openReason = enums::Lsf::NoOpenReason;
    m_openRequester = enums::Lsf::NoOpenRequest;
    m_crate = enums::Lsf::NoCrate;
    m_mode = enums::Lsf::NoMode;
    m_closeReason = enums::Lsf::NoCloseReason;
    m_closeRequester = enums::Lsf::NoCloseRequest;
    m_datagrams = 0xffffffff;
    m_modeChanges = 0xffffffff;
  }
  /// ROOT print function
  void Print(Option_t* /* option="" */) const;

  /// Create a fake for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 

  /// Compare to another in tests
  Bool_t CompareInRange( const DatagramInfo&, const std::string & name = "" ) const ; // for tests
  
private:
  
  enums::Lsf::OpenReason m_openReason;
  enums::Lsf::OpenRequester m_openRequester;
  enums::Lsf::Crate m_crate;
  enums::Lsf::Mode m_mode;    
  enums::Lsf::CloseReason m_closeReason;
  enums::Lsf::CloseRequester m_closeRequester;
  
  UInt_t m_datagrams;    
  UInt_t m_modeChanges; 

  ClassDef(DatagramInfo,1) ;

};



#endif    // ROOT_DATAGRAMINFO_H
