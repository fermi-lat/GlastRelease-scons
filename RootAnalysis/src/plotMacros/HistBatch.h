#ifndef HISTBATCH_H
#define HISTBATCH_H

#include "util.h"

class HistBatch {
private:
    static void fillListFromDir(TList *pList, TDirectory *pDir) const;
public:
    static void drawAll(const TList *l, const char *szFileName, const Int_t nCols = 1, const Int_t nRows = 1) const;
    static void drawAll(const TDirectory *d, const char *szFileName, const Int_t nCols = 1, const Int_t nRows = 1) const;
    static void drawAll(const char *szInFileName, const char *szFileName, const Int_t nCols = 1, const Int_t nRows = 1) const;
};


void HistBatch::fillListFromDir(TList *pList, TDirectory *pDir) {
    TDirectory *pOldDir = gDirectory;
    pDir->cd();
    
    TList *pKeyList = pDir->GetListOfKeys();
    
    TObject *pObj;
    TKey *pKey;

    TIter it(pKeyList);

    while (pKey = (TKey*)it()) {
        TString pClassName(pKey->GetClassName());
        if ( pClassName->BeginsWith("TH1") ||    // Cheap way of 
             pClassName->BeginsWith("TH2") ||    // not reading objects 
             pClassName->BeginsWith("TH3") )     // of types other than histos 
        {
            pObj = pKey->ReadObj();
            if (pObj->InheritsFrom("TH1")) {
                char szCycle[5];
                util::IntToStr(pKey->GetCycle(), szCycle, 5);
                ((TH1*)pObj)->SetTitle(TString(TString(((TH1*)pObj)->GetTitle()) + " [" + pDir->GetPath() + "/" + ((TH1*)pObj)->GetName() + ";" + szCycle + "]"));
                printf("\tFound histogram: %s;%d\n", ((TH1*)pObj)->GetName(), pKey->GetCycle());
                pList->Add(pObj);
            }
        } else if (pClassName == "TDirectory") {
            TDirectory *pSubDir = (TDirectory*)pKey->ReadObj();
            fillListFromDir(pList, pSubDir);
        }
    }

    pOldDir->cd();
}


void HistBatch::drawAll(const TList *pList, const char *szFileName, const Int_t nCols, const Int_t nRows) const {
    TCanvas *pCanv = new TCanvas("pCanv", "canvas", 800, 600);
    pCanv->Divide(nCols, nRows);

    TPostScript *pPs = new TPostScript(szFileName, 112);
    pPs->NewPage();

    Int_t POS_MAX = nCols * nRows, nPos = 1, nPage = 1;
    TObject *pObj;

    TIter it(pList);

    printf("\t[PAGE 1:]\n");
    while (pObj = (TObject*)it()) {
        if (nPos > POS_MAX) {
            nPos = 1;
            pCanv->Update();
            pPs->NewPage();
            printf("\t[PAGE %i]\n", ++nPage);
        }
        if (pObj->InheritsFrom("TH1")) {
            TH1* h = (TH1*)pObj;
            pCanv->cd(nPos++);
            printf("\tDrawing histogram: %s\n", ((TH1*)pObj)->GetName());
            h->Draw();
        }
    }

    pCanv->Update();
    pPs->Close();
    delete pPs;
    delete pCanv;

    printf("\tPostscript written to: %s\n", szFileName);

    // open ghostview and display file:
//    TString s("c:/gstools/gsview/gsview32.exe ");
//    s += szFileName;
//    gSystem->Exec(s);
}


void HistBatch::drawAll(const char *szInFileName, const char *szFileName, const Int_t nCols, const Int_t nRows) const {
    TDirectory *pOldDir = gDirectory;

    TFile *pFile = new TFile(szInFileName);
    TList *pList = new TList();

    fillListFromDir(pList, pFile);
    drawAll(pList, szFileName, nCols, nRows);
    
    delete pList;

    pFile->Close();
    delete pFile;
}


void HistBatch::drawAll(const TDirectory *pDir, const char *szFileName, const Int_t nCols, const Int_t nRows) const {
    TDirectory *pOldDir = gDirectory;
    pDir->cd();
    drawAll(pDir->GetList(), szFileName, nCols, nRows);
    pOldDir->cd();
}

#endif