#ifndef AcdParser_H
#define AcdParser_H 1

// $Header$

// Online EBF library includes
#include "AEMcontributionIterator.h"
#include "EBFevent.h"
// ldfReader includes
#include "ldfReader/data/TowerData.h"
// 
#include <string>

/** @class AcdParser
* @brief The specific callbacks for the CAL to handle parsing of the data.
* 
* The log method is called during the iteration of the CAL data so that the contents can
* be stored in the ldfReader object for use by ldfReader's clients.
*/

namespace ldfReader {

    class AcdParser : public AEMcontributionIterator
    {
    public:
        AcdParser(EBFevent* event,
            AEMcontribution* contribution,
            const char* prefix);

        virtual ~AcdParser() {}

        void parse();

        virtual void header(unsigned cable, AEMheader hdr);
        virtual void pha(unsigned cable, unsigned channel, ACDpha p);

        virtual int handleError(AEMcontribution *contribution, unsigned code, unsigned p1=0, unsigned p2=0) const;

        void setCurHeader(const AEMheader &header) { m_curHeader = header; };
        const AEMheader& getCurHeader() const { return m_curHeader; };

    private:
        const char *m_prefix;
        AEMheader m_curHeader;
    };
}
#endif

