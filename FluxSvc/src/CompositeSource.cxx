// $Header$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "CompositeSource.h"  //TAKE THE /.. OUT!

#include "dom/DOM_Element.hpp"
#include "facilities/Scheduler.h"
#include "facilities/SimpleEvent.h"
#include "CLHEP/Random/RandFlat.h"

// see coment below: #include "control/EventLoop.h"


#include <strstream>
#include <cassert>
#include <numeric> // for accumulate
#include <functional>
#include <iomanip>

CompositeSource::CompositeSource (double aRate)
: EventSource(aRate), m_recent(0),m_numofiters(0)//,m_time(0)
{
}

CompositeSource::~CompositeSource()
{
    for (std::vector<EventSource*>::iterator it = m_sourceList.begin();
    it != m_sourceList.end(); ++it ) delete (*it);
}


void CompositeSource::addSource (EventSource* aSource)
{
    m_sourceList.push_back(aSource);
    EventSource::setFlux( flux(EventSource::time()) );
}

void CompositeSource::rmvSource (EventSource* aSource)
{
    std::vector<EventSource*>::iterator   it = m_sourceList.begin();
    for (;it != m_sourceList.end(); ++it) {
        if ((*it) == aSource)   break;
    }
    if (it != m_sourceList.end()) {
        m_sourceList.erase(it);
        EventSource::setFlux( flux(EventSource::time()) );
    }
}

FluxSource* CompositeSource::event (double time)
{
    EventSource::setTime(time);
    
    m_numofiters=0;
    double mr = rate(EventSource::time());
    
    if( m_sourceList.size()==1 || mr ==0) {
        m_recent = m_sourceList.front();
    }else {
        
        // more than one:: choose on basis of relative rates
        double  x = RandFlat::shoot(mr), y = 0;
        std::vector<EventSource*>::iterator  now = m_sourceList.begin();
        std::vector<EventSource*>::iterator  it = now;
        
        double intrval=0.,intrmin=100000.;
        for (int q=0 ; now != m_sourceList.end(); ++now) {
            intrval=(*now)->interval(EventSource::time());
            
            if(intrval < intrmin){
                it=now;
                intrmin=intrval;
                m_numofiters=q;
            }
            //y += fabs((*it)->rate(m_time));
            //if (x <= y) {
            //    m_recent = (*it);
            //    break;
            //}
            
            m_recent = (*it);
            q++;
        }
        setInterval(intrmin);
    }
    //update the time
    //m_time += interval(m_time);
    // now ask the chosen one to generate the event, if there is a rate
    return m_recent->event(time);
}

std::string CompositeSource::fullTitle () const
{
    std::strstream  s;
    std::vector<EventSource*>::const_iterator	it = m_sourceList.begin();
    
    while (it != m_sourceList.end()) {
        
        s << (*it)->fullTitle() << " ";
        ++it;
        if (it != m_sourceList.end())    s << "+ ";
    }
    s << '\0';
    std::string t(s.str());
    s.freeze(false);
    return t;
}

std::string CompositeSource::displayTitle () const
{
    return (m_recent == 0) ? "" : m_recent->displayTitle();
}

double CompositeSource::rate(double time) const
{
    //m_time += m_time-time;
    std::vector<EventSource*>::const_iterator it = m_sourceList.begin();
    double	total_rate = 0.;
    for(;it != m_sourceList.end();++it) {
        double rr = fabs((*it)->rate(time));
        total_rate += rr;
    }
    return total_rate;
}

void	CompositeSource::setRate ( double value )
{
    double  f = rate(EventSource::time());
    if (f == 0.)    return;
    
    std::vector<float>	fvec;
    std::vector<EventSource*>::iterator it = m_sourceList.begin();
    
    while (it != m_sourceList.end())	{
        
        (*it)->setRate( value * (*it)->rate(EventSource::time())/f );
        ++it;
    }
    EventSource::setRate( value );
}

// implement virtual function
void CompositeSource::setupXML (const DOM_Element&) {}

void CompositeSource::printOn(std::ostream& out)const
{
    out << "Source(s), total rate="<< rate(EventSource::time()) << std::endl;
    
    for( std::vector<EventSource*>::const_iterator it = m_sourceList.begin();
    it != m_sourceList.end();++it)	{
        out <<  std::setw(8) << std::setprecision(4) << (*it)->rate(EventSource::time()) <<" Hz, "
            << '#' << std::setw(6) << (*it)->eventNumber() <<' '
            << (*it)->name() << ' '<< (*it)->fullTitle() << std::endl;
        
    }
    
}

std::string CompositeSource::findSource()const
{
    return m_recent->fullTitle();
}

/// return a unique number correcponding to that spectrum
int CompositeSource::numSource()const
{
    return m_numofiters;
}


/// interval to the next event
//double interval (double){
//    return m_interval;
//}

/// set the interval to the next event
//double setInterval (double interval){
//m_interval = interval;
//}


