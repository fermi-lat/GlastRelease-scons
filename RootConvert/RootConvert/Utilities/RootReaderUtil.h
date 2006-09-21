
#ifndef RootConvert_RootReaderUtil_H
#define RootConvert_RootReaderUtil_H


class StatusCode;
class TChain;

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Property.h"


namespace RootPersistence {
    
    StatusCode addFilesToChain(TChain *chain,  std::string &fileName, 
                               StringArrayProperty &fileList, MsgStream &log); 
    
    bool fileExists(std::string &fileName);
    
}

#endif
