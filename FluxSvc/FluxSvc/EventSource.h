// $Header$


#ifndef EventSource_h
#define EventSource_h 1

/** 
* \class EventSource
*
* \brief  Interface to the EventSource class. Provides a generic interface to 
* what Glastsim uses as a data source.  This can constitute GISMO based 
* event generation, reading in from a file, etc...
* 
* $Header$
*/


#include <string>

#include "../src/GPS.h"


class DOM_Element;
class FluxSource;


class EventSource
{
public:
    /// ctor/dtor
    EventSource (double aFlux = 1.0, unsigned acode = 0);
    EventSource (const DOM_Element& xelem);
    virtual ~EventSource();
    
    ///    a randomized interval to the next event - default is 1/rate()
    virtual double interval (double) = 0;
    
    ///    calculate the rate for a given flux/solid angle integral (NOTE: integral of solid angle)
    // virtual double  rate ( double solid_angle, double flux );	
    virtual double  rate (double time)const;
    virtual void    setRate ( double );
    
    ///    abstract method - create an event
    virtual FluxSource* event (double) = 0;	  
    
    ///    UI titles - used for tuple header (verbose) or window title (display)
    virtual std::string fullTitle () const;
    virtual std::string displayTitle () const;
    
    ///    flux for this source in (p/(m^2*sr*sec))
    virtual double	flux (double time) const;
    virtual void      setFlux (double value);
    
    ///    disable/enable, test this particular source 
    void      disable (){m_enabled=false;}
    void      enable(){m_enabled=true;}
    bool enabled()const{return m_enabled;}
    
    ///    integral of solid angle over which flux is incident
    virtual double	solidAngle () const;
    void      solidAngle ( double );
    
    ///	name of this flux source - for UI
    const std::string& name () const;
    void name (const std::string& value);
    
    ///    code - for monte-carlo study
    unsigned  code () const;
    virtual void  code ( unsigned );
    
    ///    area 
    static double	totalArea ();
    static void	totalArea ( double value );
    
    ///    makeOrbit - creates an orbit appropriate for this source       
    //virtual Orbit*    makeOrbit () const;
    
    /// virtual event number: should be filled in by subclass
    virtual int eventNumber()const{return -1;} 
    
    ///  say which source created the current particle
    virtual std::string findSource()const{return "";}
    
    
    /// return a unique number correcponding to that spectrum
    virtual int numSource()const{return -1;}
    
    virtual double time()const{return m_time;}
    virtual void setTime(double time){m_time=time;}
    
    
    //return how many sources are in the sourcelist (defaults to 1 if only a single FluxSource)
    virtual int howManySources(){return 1;}

    /// write the characteristics of the current source distribution to a stream
    virtual void writeSourceCharacteristic(std::ostream& out){
        out << fullTitle() << std::endl;
        out<< "default message - no sources" << std::endl;
    }
    
    //double m_time;    // elapsed time, really only needed for EventSource
private:
    double m_time;    // elapsed time, really only needed for EventSource
    
    bool m_enabled;           // toggle that it is enabled
    double m_flux;		// representative flux for this event source...
    double m_solid_angle;     // solid angle over which the rate is computed	(in steradians)
    std::string m_name;       // name of the event source (UI)
    unsigned  m_code;         // code id for this event source - for identification in the tuple.
    
    static unsigned int  s_id;    // id for new EventSources...
    static double s_total_area;   // total area for flux generation (in square meters)
};

// inline function declarations:

inline	double EventSource::solidAngle () const { return m_solid_angle; }
inline	void EventSource::solidAngle ( double value ) { m_solid_angle = value; }

inline const std::string& EventSource::name () const	{   return m_name;  }
inline void EventSource::name (const std::string& value)    { m_name = value;   }

inline double    EventSource::totalArea () { return s_total_area; }
inline void    EventSource::totalArea (double value) { s_total_area = value; }

inline unsigned EventSource::code () const { return m_code; }
inline void EventSource::code ( unsigned c ) { m_code = c; }

#endif
