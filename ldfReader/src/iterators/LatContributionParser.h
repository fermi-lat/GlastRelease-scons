#ifndef LdfReader_LatContributionParser_H
#define LdfReader_LatContributionParser_H 1

#include "LATcontributionIterator.h"
#include "EbfEventParser.h"
#include "EBFevent.h"

/** @class LatContributionIterator
@brief Provides callbacks for each component.
$Header$
*/
namespace ldfReader {
    class LatContributionParser : public LATcontributionIterator
    {
    public:
        LatContributionParser() : LATcontributionIterator()  {}
        virtual ~LatContributionParser() {}

        virtual int EBF(EBFevent* start, EBFevent* end);
    private:
        EbfEventParser m_eep;
    };

}
#endif

