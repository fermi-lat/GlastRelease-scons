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
    void drawAeff(std::string ps_filename, std::string page_title="", std::string hist_title="Effective area vs energy");

    bool fileExists();
    
    // name to use for the friend file
    std::string friend_filename();

    void open_input_file(); // override from base class

    //! the scaling factor to use for projecting scaled psf
    double psf_scale(double energy, double costh,  bool thin);
    
    static double probSum[2]; // for defining quantiles

};

#endif // RootAnalysis_irfAnalysis_PSF_h
