#ifndef RootAnalysis_irfAnalysis_PSF_h
#define RootAnalysis_irfAnalysis_PSF_h

#include "IRF.h"
#include "TProfile.h"

class PSF : public IRF {
public:
    PSF();
    void project();
    void draw(std::string ps_filename);
    void drawError(std::string ps_filename);
    void drawAsymmetry(std::string ps_filename);
    void drawAeff(std::string ps_filename);

    bool fileExists(){
        TFile f(summary_filename().c_str());
        return f.IsOpen();
    }

    static double probSum[2]; // for defining quantiles
    // histogram and display
    int nbins; 
    double xmin,xmax, ymax;

};

#endif // RootAnalysis_irfAnalysis_PSF_h
