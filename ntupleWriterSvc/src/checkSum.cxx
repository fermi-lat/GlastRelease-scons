/** 
 * @file checkSum.cxx
 *
 * @brief  checkSum is a simple class that plugs itself on RootTupleSvc to
 * calculate a checksum of the bytes being written to the MeritTuple ntuple.
 *
 * @author Michael Kuss
 *
 * $Header$
 */

#include "checkSum.h"

#include "facilities/Util.h"

#include "TLeafD.h"
#include "TLeafF.h"
#include "TLeafI.h"

#include "GaudiKernel/MsgStream.h"

#define FATAL   ( 1 )
#define ERROR   ( FATAL )
#define WARNING ( ERROR )
// e.g., to disable INFO, DEBUG, and VERBOSE:
// #define INFO    ( WARNING && 0 )
#define INFO    ( WARNING )
#define DEBUG   ( INFO )
#define VERBOSE ( DEBUG )


checkSum::checkSum(const std::string f) : m_is_open(false), m_bad(false) {
    m_fileName = f;
    facilities::Util::expandEnvVar(&m_fileName);
    if ( m_fileName.size() ) {
        m_out.open(m_fileName.c_str());
        m_is_open = m_out.is_open();
        m_bad = !m_is_open;
    }
    
}

void checkSum::write(TTree* t) {
    const TObjArray* lcol = t->GetListOfLeaves();
    const int lsize = lcol->GetEntries();
    if ( DEBUG )
        std::cout << "TTree " << t->GetName() << " has " << lsize << " leaves "
                  << std::endl;
    std::vector<unsigned char> charCol;
    Double_t EvtEventId = -1;     // initialize with something unreasonable
    Double_t EvtElapsedTime = -1;

    for ( int i=0; i<lsize; ++i ) {
        // there exists a TTreeFriendLeafIter, but how to use it?

        TObject* l = lcol->At(i);
        const std::string c = l->ClassName();
        const std::string n = l->GetName();
        if ( DEBUG )
            std::cout << i << ' ' << n << ' ' << c;
        if ( c == "TLeafD" ) {
            const Double_t v = dynamic_cast<TLeafD*>(l)->GetValue();
            if ( DEBUG )
                std::cout << ' ' << std::setprecision(25) << v
                          << std::setprecision(0) << std::endl;
            const unsigned int s = sizeof(v);
            unsigned char p[s];
            for ( unsigned int ip=0; ip<s; ++ip )
                p[ip] = 255;
            memcpy(p, &v, s);
            for ( unsigned int ip=0; ip<s; ++ip ) {
                if ( VERBOSE )
                    std::cout << "   " << i << ' ' << ip << ' '
                              << (unsigned short)p[ip] << std::endl;
                charCol.push_back(p[ip]);
            }

            if ( n == "EvtEventId" )
                EvtEventId = v;
            else if ( n == "EvtElapsedTime" )
                EvtElapsedTime = v;
        }
        else if ( c == "TLeafF" ) {
            const Float_t v = dynamic_cast<TLeafF*>(l)->GetValue();
            if ( DEBUG )
                std::cout << ' ' << std::setprecision(25) << v
                          << std::setprecision(0) << std::endl;
            const unsigned int s = sizeof(v);
            unsigned char p[s];
            for ( unsigned int ip=0; ip<s; ++ip )
                p[ip] = 255;
            memcpy(p, &v, s);
            for ( unsigned int ip=0; ip<s; ++ip ) {
                if ( VERBOSE )
                    std::cout << "   " << i << ' ' << ip << ' '
                              << (unsigned short)p[ip] << std::endl;
                charCol.push_back(p[ip]);
            }
        }
        else if ( c == "TLeafI" ) {
            const Int_t v = static_cast<Int_t>(dynamic_cast<TLeafI*>(l)->GetValue());
            if ( DEBUG )
                std::cout << ' ' << std::setprecision(25) << v
                          << std::setprecision(0) << std::endl;
            const unsigned int s = sizeof(v);
            unsigned char p[s];
            for ( unsigned int ip=0; ip<s; ++ip )
                p[ip] = 255;
            memcpy(p, &v, s);
            for ( unsigned int ip=0; ip<s; ++ip ) {
                if ( VERBOSE )
                    std::cout << "   " << i << ' ' << ip << ' '
                              << (unsigned short)p[ip] << std::endl;
                charCol.push_back(p[ip]);
            }
        }
        else {
            if ( DEBUG )
                std::cout << std::endl;
            if ( WARNING )
                std::cout << "class " << c << " is not implemented!"<<std::endl;
        }
    }
    const unsigned long theSum = simple(&charCol);

    if ( DEBUG )
        std::cout << "checksum: "
                  << std::setprecision(25)
                  << std::resetiosflags(std::ios::scientific) << EvtEventId << ' '
                  << std::setiosflags(std::ios::scientific) << EvtElapsedTime <<' '
                  << theSum << ' '
                  << charCol.size()
                  << std::endl;
    m_out.precision(25);
    m_out << std::setw(10)
          << std::resetiosflags(std::ios::scientific) << EvtEventId << "     "
          << std::setw(25)
          << std::setiosflags(std::ios::scientific) << EvtElapsedTime << "     "
          << std::setw(25) << theSum
          << std::endl;
}



unsigned long checkSum::simple(std::vector<unsigned char>* v) {
    unsigned long sum = 0;
    for( std::vector<unsigned char>::iterator it=v->begin(); it<v->end(); ++it )
        sum += *it;
    return sum;
}
