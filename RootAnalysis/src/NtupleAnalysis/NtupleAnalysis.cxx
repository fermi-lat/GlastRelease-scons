#include "NtupleAnalysis.h"
#include "../plotMacros/HistBatch.h"
#include "../plotMacros/util.h"


NtupleAnalysis::NtupleAnalysis() {
    m_nCols = 1;
    m_nRows = 1;
    m_pRootFileName = 0;
    m_pRootFile = 0;
    m_pRootTree = 0;
    m_pHistDir = 0;
    m_szPsFileName = 0;
}


NtupleAnalysis::~NtupleAnalysis() {
    recreate();
}


void NtupleAnalysis::recreate() {
    if (m_pRootFile) {
        m_pRootFile->Close();
        delete m_pRootFile;
        m_pRootFile = 0;
    }

    if (m_pRootTree) {
        m_pRootTree = 0;
    }

    if (m_pRootFileName) {
        delete m_pRootFileName;
        m_pRootFileName = 0;
    }    
}


void NtupleAnalysis::beginDraw(const char *szPsFileName) {
    m_szPsFileName = szPsFileName;
    m_szHistName[0] = 'h';
    m_nHist = 0;
    m_pHistDir = new TFile("temp.root", "RECREATE");
    m_pHistDir->cd();
}


void NtupleAnalysis::endDraw() {
    HistBatch::drawAll(m_pHistDir, m_szPsFileName, m_nCols, m_nRows);
    if (m_pHistDir) delete m_pHistDir;
}


void NtupleAnalysis::drawHist(const char *szPlotName, const TCut cut) {
    util::IntToStr(m_nHist++, m_szHistName + 1, 6);

    TString s(szPlotName);
    s += ">>";
    s += m_szHistName;
    m_pRootTree->Draw(s, cut, "");
}


Bool_t NtupleAnalysis::openRootFile(const char *szRootFileName, const char *treeName) {
    recreate();

    m_pRootFile = new TFile(szRootFileName);
    
    if (!m_pRootFile) {
        printf("Could not open file %s\n", szRootFileName);
        return kFALSE;
    }

    m_pRootTree = (TTree*)gDirectory->Get(treeName);

    if (!m_pRootTree) {
        printf("Could not open tree %s\n", treeName);
        return kFALSE;
    }

    m_pRootFileName = new TString(szRootFileName);
}


Bool_t NtupleAnalysis::drawPlots( const char *szPsFileName, const TCut cut, 
                const TString plot1,
                const TString plot2,
                const TString plot3,
                const TString plot4,
                const TString plot5,
                const TString plot6,
                const TString plot7,
                const TString plot8,
                const TString plot9,
                const TString plot10
                ) 
{
    beginDraw(szPsFileName);

    if (plot1 != "") drawHist(plot1, cut);
    if (plot2 != "") drawHist(plot2, cut);
    if (plot3 != "") drawHist(plot3, cut);
    if (plot4 != "") drawHist(plot4, cut);
    if (plot5 != "") drawHist(plot5, cut);
    if (plot6 != "") drawHist(plot6, cut);
    if (plot7 != "") drawHist(plot7, cut);
    if (plot8 != "") drawHist(plot8, cut);
    if (plot9 != "") drawHist(plot9, cut);
    if (plot10 != "") drawHist(plot10, cut);

    endDraw();
    return kTRUE;
}


//! write all plots to ps file
Bool_t NtupleAnalysis::drawAllPlots(const char *szPsFileName, const TCut cut) {
    beginDraw(szPsFileName);
    
    TObjArray *pArr = m_pRootTree->GetListOfBranches();
    for (Long64_t i = 0; i < pArr->GetEntries(); i++) {
        drawHist(((TNamed*)pArr->At(i))->GetName(), cut);
    }

    endDraw();
    return kTRUE;
}


Bool_t NtupleAnalysis::drawPlots(const char *szPsFileName, const char *szPlotFileName) {
    if (!m_pRootTree)
        return kFALSE;

    FILE *pFile = fopen(szPlotFileName, "r");
    if (!pFile)
        return kFALSE;

    TCut cut(""), cutTemp("");
    TString s = "", sTemp = "";
    char *pEnd, *pStart;
    int nCharsRead = 0, nLineLen = 0;

    int nHist = 0;
    char szHistName[7];
    szHistName[0] = 'h';

    const int BUFLEN = 255;
    char szBuffer[BUFLEN];
    char szLineBuf[BUFLEN];
    int nLineBufLen = 0;

    beginDraw(szPsFileName);

    while (nCharsRead = fread(szBuffer, sizeof(char), BUFLEN, pFile)) {
        pStart = pEnd = szBuffer;
        while (pEnd < pStart + nCharsRead) {
            if (*pEnd != '\n') {
                szLineBuf[nLineBufLen++] = *pEnd;
                pEnd++;
            } else {
				++pEnd;
                szLineBuf[nLineBufLen] = 0;
                nLineBufLen = 0;
                s = szLineBuf;
                if (s.BeginsWith("PLOT+")) {
                    s = s.Data() + 6;
                    printf("\tPLOT+ %s\n", s.Data());
                    drawHist(s.Data(), cut);
                } else if (s.BeginsWith("PLOT")) {
                    s = s.Data() + 5;
                    printf("\tPLOT %s\n", s.Data());
                    drawHist(s.Data(), cut);
                } else if (s.BeginsWith("CUT+")) {
                    printf("\tCUT+ %s\n", s.Data() + 5);
                    cut = cut || (s.Data() + 5);
                } else if (s.BeginsWith("CUT*")) {
                    printf("\tCUT* %s\n", s.Data() + 5);
                    cut = cut && (s.Data() + 5);
                } else if (s.BeginsWith("CUT")) {
                    s = s.Data() + 4;
                    printf("\tCUT %s\n", s.Data());
                    cut = s;
                } else if (s.BeginsWith("ROOT")) {
                    s = s.Data() + 5;
                    printf("\tROOT %s\n", s.Data());
                    gROOT->ProcessLineSync(s.Data());
                } else {
                    // Comment or other, do nothing
                }
            }
        }
    }

    endDraw();

    fclose(pFile);

    return kTRUE;
}


Bool_t NtupleAnalysis::makePlotFile(const char *szPlotFileName) {
    if (!m_pRootTree)
        return kFALSE;

    FILE *pFile = fopen(szPlotFileName, "w");
    if (!pFile)
        return kFALSE;

    fprintf(pFile, "# NtupleAnalysis plot file\n");
    fprintf(pFile, "# Lines have the syntax COMMAND DATA\n");
    fprintf(pFile, "# Where COMMAND is one of:\n");
    fprintf(pFile, "#   #     - comment\n");
    fprintf(pFile, "#   CUT   - replace current cut with DATA\n");
    fprintf(pFile, "#   CUT+  - (logical) OR current cut, with, DATA\n");
    fprintf(pFile, "#   CUT*  - (logical) AND current cut, with, DATA\n");
    fprintf(pFile, "#   PLOT  - plot tree branch name contained in DATA\n");
    fprintf(pFile, "#   PLOT+ - append to last plot tree branch name contained in DATA\n");
    fprintf(pFile, "#         (Plot 2D/3D histograms by separating multiple branch names with a colon)\n");
    fprintf(pFile, "#   ROOT  - pass command contained in DATA to ROOT C++ interpreter\n");
    fprintf(pFile, "#         (Useful for setting drawing styles, legend contents, etc.)\n");
    fprintf(pFile, "##################################################################################\n");
    fprintf(pFile, "# This plot file was generated from root file: %s\n", m_pRootFileName.Data());
    fprintf(pFile, "##################################################################################\n");
    fprintf(pFile, "#\n");

    TObjArray *pArr = m_pRootTree->GetListOfBranches();
    
    for (Long64_t i = 0; i < pArr->GetEntries(); i++)
        fprintf(pFile, "#PLOT %s\n", ((TNamed*)pArr->At(i))->GetName());
    
    fclose(pFile);

    return kTRUE;
}
