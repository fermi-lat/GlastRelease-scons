#ifndef RootAnalysis_irfAnalysis_PSF_h
#define RootAnalysis_irfAnalysis_PSF_h

#include "IRF.h"

class PSF : public IRF {
public:
    PSF(std::string summary_root_filename="ps.root");
    void project( double xmin, double xmax, int nbins=50);
    void draw(std::string ps_filename, double ymax, std::string title="PSF");
    void drawError(std::string ps_filename);
    void drawAsymmetry(std::string ps_filename);
    void drawAeff(std::string ps_filename);

    bool fileExists(){
        TFile f(summary_filename().c_str());
        return f.IsOpen();
    }


    std::string friend_file(){
        return output_file_root()+"psf_friend.root";
    }

    void PSF::makeFriend(); 
    
    static double probSum[2]; // for defining quantiles


};

#endif // RootAnalysis_irfAnalysis_PSF_h
