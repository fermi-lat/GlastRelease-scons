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
    //here, set the flux in a random fashion, determined by the log N/ log S characteristic graph
    double flux=1.0;
    flux = getRandomFlux();
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



double CompositeDiffuse::getRandomFlux(){
    long double prob=RandFlat::shoot(0.0, 1.0);
    long double dx=0.0000001;

    long double highthresholdofintensity = 0.0001; //for debugging and checking purposes

    long double i=0.000000001;//lowthresholdofintensity;
    
    while(prob > 0 && i<highthresholdofintensity){
        dx=0.01*pow(10.0,log10(i));
        prob-=(dx)*pofi(i);
        i+=dx;
        //	printf("\nin the findandaddnew loop; i=%12.10e, prob=%12.10e, dx=%12.10e , logi=%lf\n",i,prob,dx,log10(i));
    }
    return i;
}


long double CompositeDiffuse::pofi(long double intensity){  //this function gives P(I).  see documentation.
    long double p;
    //printf("\nabout to calculate pofi...");
    if(intensity>=pow(10.0,-7)){
        p=2.49555*pow(10.0,-13)*pow(intensity,-2.5);//egret range value for pofi
    }else{
        p=2.49555*pow(10.0,-13)*pow(intensity,-2.5-0.05*(7.0+(log(intensity)/log(10.0))));	
    }	//glast range value for pofi
    //printf("\np=%12.10e, intensity=%12.10e",p,intensity);
    return p;
}





