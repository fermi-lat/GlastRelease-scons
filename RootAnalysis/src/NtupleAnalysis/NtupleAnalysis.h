#ifndef NTUPLEANALYSIS_H
#define NTUPLEANALYSIS_H

/*! \class NtupleAnalysis
\brief Class to handle creation of PostScript files contains plots of data contained in 
a ROOT ntuple file.  
*/

class NtupleAnalysis {
private:

    TFile *m_pRootFile;
    TString *m_pRootFileName;
    TTree *m_pRootTree;

    Int_t m_nCols, m_nRows;
    Float_t m_nAspectRatio;
    Int_t m_nPsStyle;

    TDirectory *m_pHistDir;
    Int_t m_nHist;
    char m_szHistName[7];
    const char *m_szPsFileName;

    void beginDraw(const char *szPsFileName);
    void endDraw();
    void drawHist(const char *szPlotName, const TCut cut);
    void recreate();

public:

    NtupleAnalysis();
    ~NtupleAnalysis();
    
    inline void setNumCols(Int_t nCols) { m_nCols = nCols; };
    inline void setNumRows(Int_t nRows) { m_nRows = nRows; };
    inline Int_t getNumCols() { return m_nCols; };
    inline Int_t getNumRows() { return m_nRows; };

    //! select a root file for analysis
    Bool_t openRootFile(const char *szRootFileName, const char *treeName="PDR/t1");

    //! draw (up to 10) plots passed as parameters
    Bool_t drawPlots(   const char *szPsFileName, const TCut cut, 
                        const TString plot1,
                        const TString plot2 = "",
                        const TString plot3 = "",
                        const TString plot4 = "",
                        const TString plot5 = "",
                        const TString plot6 = "",
                        const TString plot7 = "",
                        const TString plot8 = "",
                        const TString plot9 = "",
                        const TString plot10 = ""
                    );
    
    //! write plots contained in plot file to ps file
    Bool_t drawPlots(const char *szPsFileName, const char *szPlotFileName = "plots.txt");


    //! write all plots to ps file, and possibily apply a cut to all data
    Bool_t drawAllPlots(const char *szPsFileName, const TCut cut="");

    //! create a plot file based on the ntuple stored in the currently loaded root file
    Bool_t makePlotFile(const char *szPlotFileName = "plots.txt");
};

#endif