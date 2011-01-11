/*
 *
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
 *
 */

#include <iomanip>

void scanTuple(const TString filename="merit.root",
               const TString treename="MeritTuple") {

    const Bool_t debugFlag = kFALSE;

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
        ULong64_t EvtEventId = -123456789;
        Float_t McId = -123456789;
        Double_t EvtElapsedTime = -123456789;

        TMD5 md5;
        for ( Int_t itl=0; itl<tlSize; ++itl ) {
            const TLeaf* tl = (TLeaf*)tlList->UncheckedAt(itl);
            const TString name = tl->GetName();
            if ( name(0,3) == "Aud" ) {
                if ( debugFlag )
                    std::cout << "Skipping auditor variable " <<name<<std::endl;
                continue;
            }
            const TString type = tl->IsA()->GetName();
            Int_t len = tl->GetLen();
            if ( len < 1 ) {
                std::cout << "Leaf " << name << " of type " << type
                          << " has length " << len << std::endl;
                continue;
            }
            // I was asking myself: "what is len?".
            // Reply: some leaves are arrays, but most of them are of
            // fundamental types and have length 1!
            for ( Int_t l=0; l<len; ++l ) {
                // TLeaf.GetValue() is of type Double_t for all leave types!
                // Thus, there is no point using the original type of the
                // contained variable, e.g. Int_t for TLeafI!
                const Double_t v = tl->GetValue(l);
                const UInt_t s = sizeof(v);
                md5.Update((UChar_t*)&v, s);
                if ( name == "EvtRun" )
                    EvtRun = v;
                else if ( name == "EvtEventId64" )
                    EvtEventId = v;
                else if ( name == "McId" )
                    McId = v;
                else if ( name == "EvtElapsedTime" )
                    EvtElapsedTime = v;
            }
            if ( debugFlag ) {
                //                std::cout << "Leaf " << name << " of type " << type
                //                          << " has value " << v << std::endl;
                std::cout << name << ' ' << type;
                if ( len != 1 )
                    std::cout << '[' << len << ']';
                for ( Int_t l=0; l<len; ++l ) {
                    std::cout << ' ';
                    tl->PrintValue(l);
                }
                std::cout << std::endl;
            }
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
