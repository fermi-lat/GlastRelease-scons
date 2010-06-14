#ifndef LATCOMPONENTPARSER_H
#define LATCOMPONENTPARSER_H 1

#include "LATcomponentIterator.h"

/** @class LatComponentParser
@brief Provides callbacks for each component.
$Header$
*/
namespace ldfReader {
    class LatComponentParser : public LATcomponentIterator
    {
    public:
        LatComponentParser() : _calSrc(0xffffffff) {}
        virtual ~LatComponentParser() {}

        virtual int UDFcomponent(EBFevent* event, EBFcontribution* contribution);
        virtual int OSWcomponent(EBFevent* event, OSWcontribution* contribution);
        virtual int GLTcomponent(EBFevent* event, GEMcontribution* contribution);
        virtual int ACDcomponent(EBFevent* event, AEMcontribution* contribution);
        virtual int CALcomponent(EBFevent* event, CALcontribution* contribution);
        virtual int TKRcomponent(EBFevent* event, TKRcontribution* contribution);
        virtual int diagnostic  (EBFevent* event, TEMcontribution* contribution);
        virtual int error       (EBFevent* event, TEMcontribution* contribution);
        virtual int cleanup     (EBFevent* event, TEMcontribution* contribution);

    private:
        unsigned _calSrc;
    };
}
#endif

