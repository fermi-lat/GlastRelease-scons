// $Header$
//  Author: T. Burnett

#ifndef __CONS_H
#define __CONS_H

#include "geometry/Cone.h"


class Cons: public Cone
{
    
public:
    Cons( double l, double ri1, double ro1, double ri2, double ro2,
        double ph1, double ph2 );
    virtual ~Cons() {}
    
    double angle1() const;
    double angle2() const;
    // parameters in addition to those of Cone
    
    virtual const char *nameOf() const { return "Cons"; }
    virtual void printOn( std::ostream& os = std::cout ) const;
    
    
};

#endif

