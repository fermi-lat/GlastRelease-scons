// $Header$
#ifndef TimeCandle_H
#define TimeCandle_H
/** 
* \class TimeCandle
*
* Spectrum: base class for energy spectrum objects
* TimeCandle: define a particle with a constant time of arrival.
* 
* $Header $
*/

#include "Spectrum.h"
#include <string>

class DOM_Element;



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//!  a convenient Spectrum : a single particle at a constant incremental time, 

class TimeCandle : public Spectrum {
public: 
    TimeCandle(const std::string& params);
    
    TimeCandle();
    //void setPosition ( float /*lat*/, float /*lon*/ ){}
    //virtual double calculate_rate(double old_rate);
    virtual float  operator()(float f)const;
    virtual const char* particleName()const;
    virtual std::string title()const;
        
    virtual std::pair<float,float> dir(float)const{
        return std::make_pair<float,float>(1.0,0.0);
    }     
    
    virtual std::pair<double,double> dir(double energy, HepRandomEngine* engine){return std::make_pair<double,double>(1.0,0.0);}
    
    double energySrc(HepRandomEngine* engine, double time){  return (*this)(engine->flat());}
    
    double interval (double time)
    {        
        return m_T0;
    }
private:
    float parseParamList(std::string input, int index);
    //float m_E0;
    double m_T0; //how many seconds pass before each next incoming particle
    std::string m_name;	// particle name to generate ("P", "gamma", ...)
    //float m_index;	// spectral index: <=1 is delta function at E0
    //float m_emax;
};

#endif
