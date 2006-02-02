#ifndef ROOT_METAEVENT_H
#define ROOT_METAEVENT_H 1

#include <iostream>
#include "TObject.h"

#include "LsfTime.h"
#include "LsfRunInfo.h"
#include "LsfDatagramInfo.h"
#include "LsfGemScalars.h"
#include "LsfConfiguration.h"

/** @class LsfMetaEvent
*
* $Header$
*/


class MetaEvent : public TObject {
  
public:
  
  MetaEvent( const RunInfo& run, const DatagramInfo& datagram, 
	     const GemScalars& scalars,
	     const Time& time,
	     const Configuration& configuration )
    :m_run(run),m_datagram(datagram),
     m_scalars(scalars),
     m_time(time),
     m_config(configuration.clone()),
     m_type(configuration.runType()){
  }
  
  MetaEvent()
    :m_config(0),
     m_type(enums::Lsf::NoRunType){
  }
  
  MetaEvent( const MetaEvent& other )
    :TObject(other),
     m_run(other.run()),
     m_datagram(other.datagram()),
     m_scalars(other.scalars()),
     m_time(other.time()),
     m_config(0),
     m_type(enums::Lsf::NoRunType){
    if ( other.configuration() != 0 ) {
      m_config = other.configuration()->clone();
      m_type = other.configuration()->runType();
    }
  }
  
  virtual ~MetaEvent(){
    delete m_config;
  }
  
  /// Information about the run this event is from
  inline const RunInfo& run() const { return m_run; }
  
  /// Information about the datagram this event came in
  inline const DatagramInfo& datagram() const { return m_datagram; }
  
  /// The extended context records
  inline const GemScalars& scalars() const { return m_scalars; }
  
  /// Information about the time markers associated with this event
  inline const Time& time() const { return m_time; } 
  
  /// Information about the configuration keys associated with this event
  inline const Configuration* configuration() const { return m_config; }
  
  /// set everything at once
  inline void initialize(const RunInfo& run, const DatagramInfo& datagram, 
			 const GemScalars& scalars,
			 const Time& time,
			 const Configuration& configuration) {
    m_run = run;
    m_datagram = datagram;
    m_scalars = scalars;
    m_time = time;
    delete m_config;
    m_config = configuration.clone();
    m_type = configuration.runType();
  }
  
  // set the individual data members
  inline void setRun( const RunInfo& val) { m_run = val; }
  inline void setDatagram( const DatagramInfo& val) { m_datagram = val; }
  inline void setScalars( const GemScalars& val) { m_scalars = val; }
  inline void setTime( const Time& val) { m_time = val; }
  inline void setConfiguration( const Configuration& configuration ) {
    delete m_config;
    m_config = configuration.clone();
    m_type = configuration.runType();
  }  

  /// Reset function
  void Clear(Option_t* /* option="" */) {
    m_run.Clear("");
    m_datagram.Clear("");
    m_scalars.Clear("");
    m_time.Clear("");
    delete m_config;
    m_config = 0;
    m_type = enums::Lsf::NoRunType;
  }

  /// ROOT print function
  void Print(Option_t* /* option="" */) const;
  
  /// Create a fake for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 
  
  /// Compare to another in tests
  Bool_t CompareInRange( const  MetaEvent&, const std::string & name = "" ) const ; 

private:
  
  /// 
  RunInfo m_run;
  DatagramInfo m_datagram;
  GemScalars m_scalars;
  Time m_time;
  Configuration* m_config;    //-> 
  
  enums::Lsf::RunType m_type;
  
  ClassDef(MetaEvent,1) ;

};

#endif    // ROOT_METAEVENT
