/** @file IBbgndTupleSelectTool.h

    @brief declaration of the IBkgndSelectTool class

$Header$

*/

#ifndef IBkgndTupleSelectTool_h
#define IBkgndTupleSelectTool_h

#include "GaudiKernel/IAlgTool.h"

#include <string>
#include <vector>

/** @class IBkgndSelectTool
    @brief manage the selection of background events to merge with signal events
    @author Dan Flath

It uses the name of a tuple variable as a key to look up trigger and downlink rates of a background
source depending on the value of that variable. See the corresponding member functions. Given such
a value, it expects to find a source of such events generated with that value which can be inserted 
into the Gleam output tuple.

*/

static const InterfaceID IID_IBkgndTupleSelectTool("IBkgndTupleSelectTool", 1 , 0);

class IBkgndTupleSelectTool : virtual public IAlgTool
{
public:

    // Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IBkgndTupleSelectTool; }

    ///! The current value of the quantity that we are selecting on
    virtual double value()const = 0;

    /** @brief select an event and copy the contents to the output tree
    */
    virtual void selectEvent() = 0;

    virtual void copyEventInfo() = 0;
 
    /**@brief the current trigger rate for the current value
    */
    virtual double triggerRate() const = 0;

    /**@brief the downlink rate for the current value
    */
    virtual double downlinkRate() const = 0;

    virtual const std::string& sourceName() const = 0;
    //virtual const std::string& name() const = 0;
};


#endif
