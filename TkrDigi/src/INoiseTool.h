/**
 * @class INoiseTool
 *
 * @brief Abstract interface to the Noise tools.
 * Currently there is but one, "General", but this for future extensions.
 *
 * @author Michael Kuss
 *
 * $Header$
 */

#ifndef __INOISETOOL_H__
#define __INOISETOOL_H__

#include "GaudiKernel/IAlgTool.h"


static const InterfaceID IID_INoiseTool("INoiseTool", 1, 0);


class INoiseTool : virtual public IAlgTool {

 public:

    /// Interface ID
    static const InterfaceID& interfaceID() { return IID_INoiseTool; }

    /**
     * Method to add noise to hits and to generate noise hits.  Returns a status
     * code upon completion.
     */
    virtual StatusCode execute() = 0;

};

#endif
