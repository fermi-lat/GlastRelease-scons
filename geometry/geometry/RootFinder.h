// $Id$

// This is a base class for any object that needs to find a root of a function

#ifndef __ROOTFINDER_H
#define __ROOTFINDER_H

class RootFinder
{
  public:

  virtual void funcd(float x, float* f, float* df)=0;
  // sub class must implement: set value of function, derivative

  float root(float x1, float x2);
  // returns the root, using the numerical recipes code for rtsafe


  static float xacc;
  static unsigned MAXIT;
  // tolerance, maximum interations
};

#endif


