#include "CalUtil/CalCalibMap.h"
#include "CalUtil/CalPedCalib.h"
#include "CalUtil/CalGainCalib.h"

template<> class CalXtalCalib<CalPedCalib,CalPedElement,2>;
template<> class CalCalibMap<CalPedCalib,CalPedElement,2>;
template<> class CalXtalCalib<CalGainCalib,CalGainElement,2>;
template<> class CalCalibMap<CalGainCalib,CalGainElement,2>;
