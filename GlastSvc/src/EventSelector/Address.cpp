// File and Version Iinformation:
// $Header$
//
// Description:
//   
//

#define Address_cpp

#include <iostream>
#include "Address.h"

Address::Address(unsigned char svc, const CLID& clid, const std::string& path)
: GenericAddress(svc, clid, path, "", 0, 0) 
{
}

