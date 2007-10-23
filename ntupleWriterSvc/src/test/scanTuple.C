/**
 * @program scanTuple
 *
 * @brief
 * small utility which scans a ntuple and calculates a checksum for each event.
 *
 * Usage: root -q -l "scanTuple.C(\"othername.root\")"
 *
 * @author Michael Kuss
 *
 * $Header$
 */

#include <iomanip>

void scanTuple(const TString filename="merit.root",
               const TString treename="MeritTuple") {

    gROOT->Reset();
    gROOT->Clear();

    //const TString filename = "merit010.root";
    //const TString treename = "MeritTuple";

    TFile* tf = new TFile(filename);
    FILE* of = fopen("checksum.txt", "w");

    if ( tf->IsZombie() ) {
        std::cerr << "Error opening file " << filename << std::endl;
        std::exit(1);
    }

    std::cout << "opening file " << filename << std::endl;
    TTree* tt = (TTree*)tf->Get("MeritTuple");
    tt->ls();
    const Int_t ttSize = tt->GetEntries();
    std::cout << "tree size: " << ttSize << std::endl;
    //    tt->SetBranchStatus("*", 1);
    //    tt->SetBranchStatus("McSourceName", 0);

    const TObjArray* tlList = tt->GetListOfLeaves();
    const Int_t tlSize = tlList->GetEntriesFast();
    std::cout << "leaflist size: " << tlSize << std::endl;

    for ( Int_t itt=0; itt<ttSize; ++itt ) {
        if( fmod(itt,100) == 0 )
            std::cout << "processed: " << (100*itt) % ttSize << std::endl;
        tt->GetEntry(itt);

        Int_t EvtRun = -123456789;
        Int_t EvtEventId = -123456789;
        Float_t McId = -123456789;
        Double_t EvtElapsedTime = -123456789;

        TMD5 md5;
        for ( Int_t itl=0; itl<tlSize; ++itl ) {
            const TLeaf* tl = (TLeaf*)tlList->UncheckedAt(itl);
            const TString name = tl->GetName();
            const TString type = tl->IsA()->GetName();
            Int_t len = tl->GetLen();
            if ( len < 1 ) {
                std::cout << "Leaf " << name << " of type " << type
                          << " has length " << len << std::endl;
                continue;
            }
            for ( Int_t l=0; l<len; ++l ) {
                if ( type == "TLeafF" ) {
                    const Float_t fv = tl->GetValue(l);
                    const UInt_t s = sizeof(fv);
                    md5.Update((UChar_t*)&fv, s);
                    if ( name == "McId" )
                        McId = fv;
                }
                else if ( type == "TLeafD" ) {
                    const Double_t dv = tl->GetValue(l);
                    const UInt_t s = sizeof(dv);
                    md5.Update((UChar_t*)&dv, s);
                    if ( name == "EvtElapsedTime" )
                        EvtElapsedTime = dv;
                }
                else if ( type == "TLeafI" ) {
                    const Int_t iv = tl->GetValue(l);
                    const UInt_t s = sizeof(iv);
                    md5.Update((UChar_t*)&iv, s);
                    if ( name == "EvtRun" )
                        EvtRun = iv;
                    if ( name == "EvtEventId" )
                        EvtEventId = iv;
                }
                else if ( type == "TLeafC" ) {
                    // GetValueString() makes the ENTIRE code very slow!
                    //const TString sv = ((TLeafC*)tl)->GetValueString();
                    const char* sv = tl->GetValuePointer();
                    const UInt_t s = len;
                    md5.Update((UChar_t*)sv, s);
                    len = 1;
                    break;
                }
                else {
                    std::cout << "Leaf type " << type
                              << " not implemented yet!" << std::endl;
                }
            }
            /*
              std::cout << type << "[" << len << "] " << name;
              for ( Int_t l=0; l<len; ++l ) {
              std::cout << ' ';
              tl->PrintValue(l);
              }
              std::cout << std::endl;
            */
        }

        md5.Final();
        printf("%7d %7d %7.0f %23.17e ", EvtRun,EvtEventId,McId,EvtElapsedTime);
        md5.Print();
        fprintf(of, "%7d %7d %7.0f %23.17e %s\n", EvtRun, EvtEventId, McId,
                EvtElapsedTime, md5.AsString());
    }

    tf->Close();
    fclose(of);
}
