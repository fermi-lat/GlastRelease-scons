#ifndef _GlastSvc_RandomAccess_H
#define _GlastSvc_RandomAccess_H 1


// Include files
#include "GlastSvc/GlastRandomSvc/IRandomAccess.h"

#include "GaudiKernel/AlgTool.h"

#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RandGauss.h"

/** @class RandomAccess
*
* @brief Gaudi Tool to set the address of CLHEP random engines
* 
* This tool contains methods to set the CLHEP random engine 
* associated with GLAST shared libraries that use CLHEP random 
* number generators It is meant to
* be inherited by a tool (generally containing no additional methods)
* created in packages that generate shared libraries, one or more
* of whose components use CLHEP random numbers. GlastRandomSvc, in
* its initialize method, queries for the existence of tools that 
* inherit from IRandomAccess and then gives each a different instance 
* of the CLHEP engine  
*     
*
* @authors Toby Burnett, Karl Young
*
* $Header$
*/

class RandomAccess : public AlgTool, virtual public IRandomAccess 
{
public:

    RandomAccess( const std::string& type, const std::string& name,
        const IInterface* parent) : AlgTool( type, name, parent ) 
    {
        // Declare additional interface
        declareInterface<IRandomAccess>(this);
    }

    //! Pass on a pointer to a new engine to be used by the local static instance of HepRandom
    //! return pointer to the old one 
    virtual HepRandomEngine* setTheEngine( HepRandomEngine* newEngine){
        HepRandomEngine* oldEngine = HepRandom::getTheEngine();
        HepRandom::setTheEngine(newEngine);
        return oldEngine;
    }

    typedef void(*SetFlag)(bool);
    virtual  SetFlag getRandSet(){
        return &RandGauss::setFlag;
    }

};

#endif // _GlastSvc_RandomAccess_H
