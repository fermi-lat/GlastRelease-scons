/**
 * @class UBinterpolateTool
 *
 * @brief interpolation class for the unbiased energy algorithm
 *
 * @author Michael Kuss, Carmelo Sgro'
 *
 * 
 */

#ifndef __UBINTERPOLATETOOL_H__
#define __UBINTERPOLATETOOL_H__

#include "GaudiKernel/IAlgTool.h"



#include <string>


// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IUBinterpolateTool("IUBinterpolateTool", 0 , 0); 


class IUBinterpolateTool : virtual public IAlgTool

{
 public:

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IUBinterpolateTool; }
    /// load map
    virtual void addBiasMap(std::string mapName, std::string calibFileName) = 0;
    /// returns an interpolated value for zDir and LogE
    virtual  float interpolate( std::string mapName,  float logE,  float zDir)  = 0;


};

#endif
