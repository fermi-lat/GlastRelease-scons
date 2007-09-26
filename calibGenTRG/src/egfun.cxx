#include "egfun.h"
#include <math.h>

Double_t egfun(Double_t *t, Double_t *parm){
  return parm[0] * (1 + (t[0] - parm[1])/parm[2]) * exp(-(t[0]-parm[1])/parm[2]);
    }
