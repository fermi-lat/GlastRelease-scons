

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"

#ifdef WIN32
#include "TSystem.h" // To get TreePlayer loaded
#endif

#include "RootConvert/Utilities/RootReaderUtil.h"

#if ROOT_VERSION(5,14,0) <=  ROOT_VERSION_CODE
#include "TStreamerElement.h"
void FixElement(TStreamerElement *el) {
   if (el) {
      el->SetTypeName("vector<double>");
      el->Update(TClass::GetClass("vector<Double32_t>"),gROOT->GetClass("vector<double>"));
   }
}
#endif

namespace RootPersistence {

    // add root files to a chain

    StatusCode addFilesToChain(TChain *chain, std::string &fileName, 
                               StringArrayProperty &fileList, MsgStream &log, bool recFile)
    {
        
        std::string emptyStr("");
        
        if (fileName.compare(emptyStr) != 0) {
            if ( !fileExists(fileName, recFile)) {
                log << MSG::ERROR << "ROOT file " << fileName.c_str()
                    << " could not be opened for reading." << endreq;
                return StatusCode::FAILURE;
            }
            chain->Add(fileName.c_str());
            log << MSG::INFO << "Opened file: " << fileName.c_str() << endreq;
        } else {
            const std::vector<std::string> fList = fileList.value( );
            std::vector<std::string>::const_iterator it;
            std::vector<std::string>::const_iterator itend = fList.end( );
            for (it = fList.begin(); it != itend; it++) {
                std::string theFile = (*it);
                if ( !fileExists(theFile, recFile)) {
                    log << MSG::ERROR << "ROOT file " << theFile.c_str()
                        << " could not be opened for reading." << endreq;
                    return StatusCode::FAILURE;
                }
                chain->Add(theFile.c_str());
                log << MSG::INFO << "Opened file: " << theFile.c_str() << endreq;
            }
        }
        
        return StatusCode::SUCCESS;
    }
            
      
    // Check if a root file can be opened.
    
    bool fileExists(std::string &filename, bool recFile)
    {
        TFile *f = TFile::Open(filename.c_str());
        if (!f) { 
            return false;
        } 
        if (f->IsZombie()) {
            delete f;
            return false;
        } else {
            if (recFile) {
#if ROOT_VERSION(5,14,0) <=  ROOT_VERSION_CODE
                if (f->GetVersion()<51200) {
                    TClass *cl = TClass::GetClass("AcdRecon");
                    FixElement((TStreamerElement*)cl->GetStreamerInfo(14)->GetElements()->FindObject("m_rowDocaCol"));
                    FixElement((TStreamerElement*)cl->GetStreamerInfo(14)->GetElements()->FindObject("m_rowActDistCol"));
                    FixElement((TStreamerElement*)cl->GetStreamerInfo(14)->GetElements()->FindObject("m_rowActDistCol_down"));
                    FixElement((TStreamerElement*)cl->GetStreamerInfo(14)->GetElements()->FindObject("m_energyCol"));
                    cl->GetStreamerInfo(14)->Clear("build");
                }
#endif
            }

            f->Close();
            delete f;
            return true;
        }
            
    }   
    
}

