// Spectrum.cxx: implementation of the Spectrum class.
//
//////////////////////////////////////////////////////////////////////

#include "FluxSvc/Spectrum.h"
#include <cmath>

// CLHEP
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Spectrum::~Spectrum(){}

//float Spectrum::fraction(float e)
//{
//    // this is essentially rtbis of numerical recipies
//    static int JMAX=40;
//    static double xacc = e*1e-3;
//    static double zero = 1e-4;

//    double dx,xmid, x1=0,x2=1;
//    double f = (*this)(x1)-e;
//    double fmid = (*this)(x2)-e;

//    // Check to see if the x2 contains the root
//    if(abs(fmid) < zero) return x2;
//
//    if (f*fmid >0.0) return 0.5; // fail?
//
//    double rtb = f < 0.0 ? (dx=x2-x1, x1) : (dx = x1-x2,x2);
//
//    for(int j=0; j< JMAX; j++) {
//	fmid = (*this)(xmid=rtb+ (dx*=0.5))-e;
//	if( fmid <0.0 ) rtb = xmid;
//	if( abs(dx) < xacc || fmid ==0) return rtb;
//    }
//    return rtb; // fail
//}

double Spectrum::flux (double time ) const {
  return 0.; // flag that we don't have a flux
}

double Spectrum::solidAngle( )const
{
    return 1.0; // flag that doesn't calculate
}

 std::pair<float,float> Spectrum::dir(float energy)const
{

    // return solid angle pair (costh, phi) for the given energy
    // default: random except for Earth occultation
    //return std::make_pair(static_cast<float>(RandFlat::shoot(-0.4,1.0)),
        //static_cast<float>(RandFlat::shoot(0.0, 2*M_PI)) );
	 //here's an attempt at random default distributions as above:
	 return std::make_pair(((rand()/32767.0)*1.4)-0.4,(rand()/32767.0)*2*M_PI);
    
}



 const char * Spectrum::particleName()const{
	 static const char x='p';
	 return &x;
 }

double Spectrum::energySrc(HepRandomEngine* engine)
{
   // default implementation, which works for other Spectrum objects
	return (*this)(engine->flat());
}

std::pair<double,double> Spectrum::dir(double energy, HepRandomEngine* engine)
{
	// default that needs fixing!
	return dir(energy);
}

void Spectrum::parseParamList(std::string input, std::vector<float>& output) const
{

	int i=0;
	for(;!input.empty() && i!=std::string::npos;){
		float f = ::atof( input.c_str() );
		output.push_back(f);
		i=input.find_first_of(",");
		input= input.substr(i+1);
	} 
}