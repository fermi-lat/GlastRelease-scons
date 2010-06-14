#ifndef ldfReader_LATCONTRIBUTIONPARSER_CXX
#define ldfReader_LATCONTRIBUTIONPARSER_CXX 1

/** @file LatContributionParser.cxx
@brief Implementation of the LatContributionParser class

$Header$
*/

//#include <stdio.h> // included for LATcomponentIterator.h in Online/EBF
#include "LatContributionParser.h"
#include "../EbfDebug.h"

namespace ldfReader {
    int LatContributionParser::EBF(EBFevent* event, EBFevent* end) {
        // Iterate over a list of EBF events
        m_eep.iterate(event, end);

        return m_eep.status();
    }

}
#endif
