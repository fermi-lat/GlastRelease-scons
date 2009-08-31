#ifndef AdfParser_H
#define AdfParser_H 1

// $Header$

/** @class AdfParser
* @brief The parsing of the Ancillary BeamTest contributions.
*
*/

namespace ldfReader {

    class AdfParser
    {
    public:
        AdfParser() { };

        virtual ~AdfParser() {}

        int parseHeader(const unsigned char* buf);
        int parseEvent(const unsigned char* buf);
        int parseTrailer(const unsigned char* buf);

    //virtual int handleError(CALcontribution*, unsigned code, unsigned p1=0, unsigned p2=0) const;

    private:
    };
}
#endif

