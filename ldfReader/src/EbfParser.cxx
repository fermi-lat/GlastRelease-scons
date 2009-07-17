#ifndef ldfReader_EbfParser_CXX
#define ldfReader_EbfParser_CXX

#include "EbfDebug.h"

/** @file EbfParser.cxx
@brief Implementation of the EbfParser base class

$Header$
*/

#include "ldfReader/data/LatData.h"
#include "ldfReader/EbfParser.h"
namespace ldfReader {

int EbfParser::setDebug(int on) {
    return EbfDebug::setDebug(on);
}

int EbfParser::setAcdRemap(const std::string& filename) {

    if (filename != "") {
       return (ldfReader::LatData::instance()->setAcdRemap(filename));
    }
    return -1;
}


void EbfParser::setIgnoreSegFault() {
   ldfReader::LatData::instance()->setIgnoreSegFault(true);

}

void EbfParser::setOldStyleRunId(bool flag) {
   ldfReader::LatData::instance()->setOldStyleRunId(flag);
}

} //namespace
#endif
