///
///    GRBUniverse: access 2 some cosmological quantities 
///    Authors: Nicola Omodei & Johann Coen Tanugi 
///
#include <math.h>
#include "GRBConstants.h"
#include "GRBUniverse.h"

GRBUniverse::GRBUniverse() {
  
  double _wzel=0.0;
  _qo=(1.0+3.0*_wzel)/2.0;
  _Hubble=cst.Hubble*1.e+5;           // cm/sec/Mpc
  _Dist=(cst.c/_Hubble/pow(_qo,2.0))*(cst.red*_qo+(_qo-1.0)*(-1.0+sqrt(2.0*_qo*cst.red+1.0)))*cst.mpc2cm;
  _Area=(4.*cst.pi)*pow(_Dist,2); // [cm^2]
}
