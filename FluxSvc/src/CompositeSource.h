// $Header$

#ifndef CompositeSource_h
#define CompositeSource_h 1

#include "FluxSvc/EventSource.h"
#include <vector>

class FluxSource;

//! holds multiple Eventsource objects ; acts as a container for them.
class CompositeSource : public EventSource { 
public:
    ///    constructor/destructor
    CompositeSource (double aRate = 1.0);
    virtual ~CompositeSource();
    
    
    ///    add a source to the list
    void addSource (EventSource* aSource);
    void rmvSource (EventSource* aSource);
    
    /// generate an event from from one of the sources 
    /// which make up the composite, and return a pointer to it
    virtual FluxSource* event (double time);
    
    /// rate - compute overall rate...
    virtual double rate (double time)const;
    virtual void   setRate ( double );
    
    /// flux into 1 m^2 integrated over angles
    virtual double flux(double time)const{
        return rate(time)/totalArea();}
    
    ///    full-length title description of this EventSource.
    virtual std::string fullTitle () const;
    
    ///    brief title description (for display) for this event source.
    virtual std::string displayTitle () const;
    
    /// dump current list of sources, rates
    void printOn(std::ostream& out)const;
    
    ///  say which source created the current particle
    std::string findSource()const;
    
    /// return a unique number correcponding to that spectrum
    int numSource()const;
    
    
    ///	    list of sources which make up this composite
    std::vector< EventSource* >& sourceList ();
    const std::vector< EventSource* >& sourceList () const;
    void sourceList (const std::vector< EventSource* >& value);
    
    /// interval to the next event
    double interval (double){return m_interval;}
    
    /// set the interval to the next event
    double setInterval (double interval){return (m_interval = interval);}
    
    // double m_time; 
    double m_interval;
    
    //return how many sources are in the sourcelist
    int howManySources(){return m_sourceList.size();}
    
protected:
    virtual void setupXML (const DOM_Element&);
    
    //number of times we've iterated the front() pointer into sourcelist 
    //to get the current particle - represents the source
    int m_numofiters;
    
    
    //private: 
    std::vector< EventSource* > m_sourceList;
    EventSource*  m_recent;
};

inline std::vector< EventSource* >& CompositeSource::sourceList ()
{
    return m_sourceList;
}

inline const std::vector< EventSource* >& CompositeSource::sourceList () const
{
    return m_sourceList;
}

inline void CompositeSource::sourceList (const std::vector< EventSource* >& value)
{
    m_sourceList = value;
}
#endif
