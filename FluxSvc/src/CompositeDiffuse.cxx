#include "CompositeDiffuse.h"
#include "CLHEP/Random/RandFlat.h"
#include "FluxSvc/FluxSource.h"
#include "SimpleSpectrum.h"

void CompositeDiffuse::addSource (EventSource* aSource)
{
m_sourceList.push_back(aSource);
double flux = aSource->flux(EventSource::time());
//EventSource::setFlux( flux );
m_unclaimedFlux-=flux;
}
/*
  void CompositeDiffuse::rmvSource (EventSource* aSource)
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
*/

FluxSource* CompositeDiffuse::event (double time)
{
    EventSource::setTime(time);
    
    int m_numofiters=0;

    // here we should be setting the total rate as the maximum sum rate of everything - FIX!
    //double mr = rate(EventSource::time());
    double mr = m_totalFlux;
    
    // do this once if there is just one source, or no rate at all (null source?)
    if( m_sourceList.size()==1 || mr == 0) {
        m_recent = m_sourceList.front();
    }else {
        
        // more than one:: choose on basis of relative rates
        double  x = RandFlat::shoot(mr), y = 0;
        std::vector<EventSource*>::iterator  now = m_sourceList.begin();
        std::vector<EventSource*>::iterator  it = now;
        
        double intrval=0.,intrmin=100000.;
        
        for ( int q=0; now != m_sourceList.end(); ++now) {
            (*now)->event(time); // to initialize particles, so that the real interval for the particle is gotten.
            intrval=(*now)->interval(EventSource::time());
            
            if(intrval < intrmin){
                it=now;
                intrmin=intrval;
                m_numofiters=q;
            }
            
            m_recent = (*it);
            q++;
        }
        //so now m_recent is the "soonest" source, and intrmin is the "soonest": time.
        // but, what if the "leftover" flux sends a source even faster?
        intrval=remainingFluxInterval();
        if(intrval < intrmin){
            
            addNewSource();
            intrmin=intrval;
            //++now;
            //m_recent = *(m_sourceList.end());
            it = m_sourceList.end();
            --it;
            m_recent = (*it);
        }
        //set the interval that won out after everything.
        setInterval(intrmin);
    }
    //update the time
    //m_time += interval(m_time);
    // now ask the chosen one to generate the event, if there is a rate
    return (FluxSource*)m_recent;//->event(time);
}


void CompositeDiffuse::addNewSource(){
    double flux=1.0;
    //Set up the spectrum with the appropriate power spectrum
    Spectrum * spec=new SimpleSpectrum("gamma", 100.0);
    
    // Now make the random vector correesponding to random galactic direction
    double  costh = -RandFlat::shoot(-1.,1./*_minCos, _maxCos*/);
    double  sinth = sqrt(1.-costh*costh);
    double  phi = RandFlat::shoot(0.0, 2*M_PI);
    
    Vector launchDir = Vector(cos(phi)*sinth, sin(phi)*sinth, costh);
    
    //Then call the FluxSource Constructor to use the galactic direction specified by launchDir.
    EventSource* aSource=new FluxSource(flux,spec,&launchDir);
    //then add it into the list of sources..
    addSource(aSource);
    //..and subtract the total flux from what remains...
    m_unclaimedFlux -= flux;
}


double CompositeDiffuse::remainingFluxInterval(){
    
    double  r = (solidAngle()*(m_unclaimedFlux)*6.);
    
    if (r <= 0){ return 1E20;
    }else{  
        double p = RandFlat::shoot(1.);
        return (-1.)*(log(1.-p))/r;
    }
    
}