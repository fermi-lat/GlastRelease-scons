// $Header$

// Implement the proton cosmic ray spectrum function


// local namespace
namespace {
    
    // define the fromMeV function which will convert from
    // MeV to GeV
    float   fromMeV ( float e ) { return e / 1000.; }
}

#include "FluxSvc/ProtonSpectrum.h"

#include <algorithm>

float ProtonSpectrum::operator()(float r)const
{   return fromMeV( protonenergy( std::max(1e-6, 1.0-r)) );
}


#include <cmath>
#include <algorithm>
#include <functional>

#include "CLHEP/Random/Random.h"
#include "FluxSvc/GPS.h"


#include "FluxSvc/SpectrumFactory.h"

static SpectrumFactory<ProtonSpectrum> factory;
const ISpectrumFactory& ProtonSpectrumFactory = factory;

//==============================================================================
// following slightly modified by THB from protonenergy.c:
//   * convert to C++, using static, declarations when appropriate
//   * pass random float [0,1] to it.
// note that it apparently returns MeV.
//------------------------------------------------------------------------------
//	     no solar flares included
//
//   Author: Patrick Nolan, Stanford University
//	September 1993
//
//   Input: none.
//
//  Output: The returned value is the energy of a proton randomly
//       chosen from the predicted spectrum in orbit.  The value
//       will be different for each call.

struct  tbl {float flux; float energy;};

int ProtonSpectrum::prolocate(struct tbl *data, int n, double x)
/*    specialized implementation of simple binary search
algorithm taken from Numerical Recipes
*/
{
    int ju, jm, jl;
    int ascnd;
    
    jl=0; ju=n+1;
    ascnd = (data[n-1].flux > data[0].flux);
    while (ju-jl > 1) {
        jm=(ju+jl) >>1;
        if ((x > data[jm-1].flux) == ascnd)
            jl = jm;
        else
            ju = jm;
    }
    return jl-1;
}

#ifdef _MSC_VER
#  pragma warning(disable:4305) //warning C4305: 'initializing' : truncation from 'const double' to 'float'
#endif //_MSC_VER

double ProtonSpectrum::protonenergy(float r)
{
    static struct  tbl
        datain[] = {
        /* integral flux, energy pairs from program  SPEC */
            { 0.78953E+01 , 0.87096E+05 },
            { 0.87333E+01 , 0.81283E+05 },
            { 0.96432E+01 , 0.75857E+05 },
            { 0.10631E+02 , 0.70794E+05 },
            { 0.11704E+02 , 0.66069E+05 },
            { 0.12869E+02 , 0.61659E+05 },
            { 0.14133E+02 , 0.57544E+05 },
            { 0.15505E+02 , 0.53703E+05 },
            { 0.16993E+02 , 0.50118E+05 },
            { 0.18609E+02 , 0.46773E+05 },
            { 0.20361E+02 , 0.43651E+05 },
            { 0.22262E+02 , 0.40738E+05 },
            { 0.24323E+02 , 0.38019E+05 },
            { 0.26557E+02 , 0.35481E+05 },
            { 0.28980E+02 , 0.33113E+05 },
            { 0.31604E+02 , 0.30903E+05 },
            { 0.34447E+02 , 0.28840E+05 },
            { 0.37526E+02 , 0.26915E+05 },
            { 0.40858E+02 , 0.25119E+05 },
            { 0.44465E+02 , 0.23442E+05 },
            { 0.48365E+02 , 0.21878E+05 },
            { 0.52583E+02 , 0.20417E+05 },
            { 0.57140E+02 , 0.19055E+05 },
            { 0.62062E+02 , 0.17783E+05 },
            { 0.67375E+02 , 0.16596E+05 },
            { 0.73107E+02 , 0.15488E+05 },
            { 0.79286E+02 , 0.14454E+05 },
            { 0.85763E+02 , 0.13490E+05 },
            { 0.92247E+02 , 0.12589E+05 },
            { 0.98383E+02 , 0.11749E+05 },
            { 0.10390E+03 , 0.10965E+05 },
            { 0.10863E+03 , 0.10233E+05 },
            { 0.11257E+03 , 0.95499E+04 },
            { 0.11599E+03 , 0.89125E+04 },
            { 0.11906E+03 , 0.83176E+04 },
            { 0.12183E+03 , 0.77624E+04 },
            { 0.12435E+03 , 0.72443E+04 },
            { 0.12666E+03 , 0.67608E+04 },
            { 0.12878E+03 , 0.63096E+04 },
            { 0.13073E+03 , 0.58884E+04 },
            { 0.13247E+03 , 0.54954E+04 },
            { 0.13394E+03 , 0.51286E+04 },
            { 0.13507E+03 , 0.47863E+04 },
            { 0.13585E+03 , 0.44668E+04 },
            { 0.13635E+03 , 0.41687E+04 },
            { 0.13668E+03 , 0.38904E+04 },
            { 0.13691E+03 , 0.36308E+04 },
            { 0.13704E+03 , 0.33884E+04 },
            { 0.13710E+03 , 0.31623E+04 },
            { 0.13712E+03 , 0.29512E+04 },
    };
    
    static int npts = sizeof(datain)/sizeof(datain[0]);  /* size of array */
    double f = r*datain[npts-1].flux;
    if (f < datain[0].flux )
        return datain[0].energy *
        pow(f/ datain[0].flux, -0.61);
    else {
        int n = prolocate(datain, npts, f);
        return datain[n].energy + (datain[n+1].energy-datain[n].energy) *
            (f-datain[n].flux) / (datain[n+1].flux-datain[n].flux);
    }
}


double ProtonSpectrum::calculate_rate(double old_rate)
{
    return old_rate;
}

const char* ProtonSpectrum::particleName()const {
    return "p";
}

std::string ProtonSpectrum::title()const {
    return "ProtonSpectrum";
}


