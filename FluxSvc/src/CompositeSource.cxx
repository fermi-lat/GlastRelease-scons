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
    //here, set up the associated vectors by default.
    m_unusedSource.push_back(0);
    m_sourceInterval.push_back(-1.);
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
    int i=0; //for iterating through the m_unusedSource vector
    int winningsourcenum; //the number of the "winning" source
    
    EventSource::setTime(time);
    
    m_numofiters=0;
    double mr = rate(EventSource::time());
    
    if( m_sourceList.size()==1 || mr ==0) {
        m_recent = m_sourceList.front();
    }else {
        
        // more than one:: choose on basis of relative rates
        // NOT used? THB  double  x = RandFlat::shoot(mr), y = 0;
        std::vector<EventSource*>::iterator  now = m_sourceList.begin();
        std::vector<EventSource*>::iterator  it = now;
        
        double intrval=0.,intrmin=100000.;
	int q;
        for (q=0 ; now != m_sourceList.end(); ++now) {
            if(m_unusedSource[i]==1){
                intrval=m_sourceInterval[i];
                //std::cout << i << " is unused, interval is "<< intrval << std::endl;
            }else{
                (*now)->event(time); // to initialize particles, so that the real interval for the particle is gotten.
                intrval=(*now)->interval(EventSource::time()); //this picks out the interval of each source
                m_unusedSource[i]=1;
                m_sourceInterval[i]=intrval;
            }
            
            if(intrval < intrmin){
                //the present source is "winning" here
                it=now;
                intrmin=intrval;
                m_numofiters=q;
                winningsourcenum=i;
            }
            
            m_recent = (*it);
            q++;
            i++;
        }
        setInterval(intrmin);
        now = m_sourceList.begin();
        for (q=0 ; now != m_sourceList.end(); ++now) {
            //this loop sets the intervals back in accordance with
            //how far ahead time will move.
            //std::cout << "decrementing source " << q << " ,  by " << intrmin << std::endl;
            m_sourceInterval[q] = m_sourceInterval[q] - intrmin;
            q++;
        }
    }
    m_unusedSource[winningsourcenum]=0; //the current "winning" source is getting used..
    // now ask the chosen one to return the event.
    return (FluxSource*)m_recent;
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

int CompositeSource::numSource()const
{
    ///Purpose: Return a unique number correcponding to the current spectrum.
    return m_numofiters;
}

