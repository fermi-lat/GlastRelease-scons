// $Header$
#ifndef LHCBEVENT_RUN_H
#define LHCBEVENT_RUN_H 1


//Include files
#include <iostream>
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/DataObject.h"
#include "Event/Utilities/TriggerPattern.h"
#include "Event/Utilities/RandomNumberSeed.h"
#include "Event/TopLevel/Definitions.h"


// Externals 
extern const CLID& CLID_Run;


//------------------------------------------------------------------------------
//
// ClassName:   Run
//  
// Description: Essential information of the run
//
//              It conttains
//              - run number
//              - run type
//              - trigger type
//              - enabled trigger mask
//              - enabled detector mask
//              - luminocity
//              - fill number
//              - generator type
//              - random number seed
//
// Author:      Pavel Binko
// Changes:     P.Binko 19/10/1999 : Formating of ASCII output
//
//------------------------------------------------------------------------------


class Run : public DataObject                                                  {

public:
  /// Constructors
  Run(const char* name = "Event")
    : DataObject(name),
      m_runType(0),
      m_triggerType(0),
      m_enabledDetectorMask(0),
      m_luminosity(0.),
      m_fillNumber(0),
      m_generatorType(0)                                                     { }
  /// Destructor
  virtual ~Run()                                                             { }

  /// Retrieve reference to class definition structure
  virtual const CLID& clID() const                    { return Run::classID(); }
  static const CLID& classID()                              { return CLID_Run; }

  /// Retrieve run type
  long runNumber () const                                                      {
    return m_run;
  }
  /// Update run type
  void setRunNumber (long value)                                               {
    m_run = value;
  }

  /// Retrieve run type
  long runType () const                                                        {
    return m_runType;
  }
  /// Update run type
  void setRunType (long value)                                                 {
    m_runType = value;
  }

  /// Retrieve trigger type
  long triggerType () const                                                    {
    return m_triggerType;
  }
  /// Update trigger type
  void setTriggerType (long value)                                             {
    m_triggerType = value;
  }

  /// Retrieve trigger enable mask
  const TriggerPattern& enabledTriggerPattern () const                         {
    return m_enabledTriggerMask;
  }
  /// Update trigger enable mask
  void setEnabledTriggerPattern (const TriggerPattern& value)                  {
    m_enabledTriggerMask=value;
  }

  /// Retrieve detector mask
  long enabledDetectorMask () const                                            {
    return m_enabledDetectorMask;
  }
  /// Update detector mask
  void setEnabledDetectorMask (long value)                                     {
    m_enabledDetectorMask = value;
  }

  /// Retrieve integrated luminosity
  double luminosity () const                                                   {
    return m_luminosity;
  }
  /// Update integrated luminosity
  void setLuminosity (double value)                                            {
    m_luminosity = value;
  }

  /// Retrieve fill number
  long fillNumber () const                                                     {
    return m_fillNumber;
  }
  /// Update fill number
  void setFillNumber (long value)                                              {
    m_fillNumber = value;
  }

  /// Retrieve MC generator type
  long generatorType () const                                                  {
    return m_generatorType;
  }
  /// Update MC generator type
  void setGeneratorType (long value)                                           {
    m_generatorType = value;
  }

  /// Retrieve initial MC random number seed
  const RandomNumberSeed randomNumberSeed () const                             {
    return m_randomNumberSeed;
  }
  /// Update initial MC random number seed
  void setRandomNumberSeed (RandomNumberSeed value)                            {
    m_randomNumberSeed = value;
  }

  /// Serialize the object for writing
  virtual StreamBuffer& serialize( StreamBuffer& s ) const;
  /// Serialize the object for reading
  virtual StreamBuffer& serialize( StreamBuffer& s );

  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s, const Run& obj )          {
    return obj.fillStream(s);
  }
  /// Fill the output stream (ASCII)
  virtual std::ostream& fillStream( std::ostream& s ) const;

private:
  /// Run number
  long                 m_run;
  /// Run type
  long                 m_runType;
  /// Trigger type
  long                 m_triggerType;
  /// Trigger enable mask
  TriggerPattern       m_enabledTriggerMask;
  /// Detector enable mask
  long                 m_enabledDetectorMask;
  /// Integrated luminosity
  double               m_luminosity;
  /// Fill number
  long                 m_fillNumber;
  /// MC generator type
  long                 m_generatorType;
  /// Initial MC random number seed
  RandomNumberSeed     m_randomNumberSeed;

};


//
// Inline code must be outside the class definition
//


/// Serialize the object for writing
inline StreamBuffer& Run::serialize( StreamBuffer& s ) const                   {
  DataObject::serialize(s);
  return s
    << m_run
    << m_runType
    << m_triggerType
    << m_enabledTriggerMask
    << m_enabledDetectorMask
    << m_luminosity
    << m_fillNumber
    << m_generatorType
    << m_randomNumberSeed;
}


/// Serialize the object for reading
inline StreamBuffer& Run::serialize( StreamBuffer& s )                         {
  DataObject::serialize(s);
  return s
    >> m_run
    >> m_runType
    >> m_triggerType
    >> m_enabledTriggerMask
    >> m_enabledDetectorMask
    >> m_luminosity
    >> m_fillNumber
    >> m_generatorType
    >> m_randomNumberSeed;
}


/// Fill the output stream (ASCII)
inline std::ostream& Run::fillStream( std::ostream& s ) const                  {
  return s
    << "class Run :\n"
    << "\n    Run number           = "
    << GlastEventField( GlastEvent::field12 )
    << m_run
    << "\n    Run type             = "
    << GlastEventField( GlastEvent::field12 )
    << m_runType
    << "\n    Trigger type         = "
    << GlastEventField( GlastEvent::field12 )
    << m_triggerType
    << "\n    Present triggers     = " << m_enabledTriggerMask
    << "\n    Detector enable mask = "
    << GlastEventField( GlastEvent::field12 )
    << m_enabledDetectorMask
    << "\n    Luminosity           = "
    << EventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_luminosity
    << "\n    Fill number          = "
    << GlastEventField( GlastEvent::field12 )
    << m_fillNumber
    << "\n    MC generator type    = "
    << GlastEventField( GlastEvent::field12 )
    << m_generatorType
    << "\n    Random number seed   = " << m_randomNumberSeed;
}


#endif // LHCBEVENT_RUN_H
