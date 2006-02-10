#ifndef ldfReader_EbfParser_CXX
#define ldfReader_EbfParser_CXX

#include "EbfDebug.h"

/** @file EbfParser.cxx
@brief Implementation of the EbfParser base class

$Header$
*/

#include "ldfReader/EbfParser.h"
namespace ldfReader {

bool EbfParser::setDebug(bool on) {
    return EbfDebug::setDebug(on);
}

}


#endif
