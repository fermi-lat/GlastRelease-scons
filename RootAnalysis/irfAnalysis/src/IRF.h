/** @file IRF.h
  @brief Create a set of histograms to allow analysis of the  response

  $Header$
*/
#ifndef IRF_H
#define IRF_H
#include "TCut.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFrame.h"

#include <stdio.h>
#include <string>
#include <iostream>

/**
    @class IRF
    @brief Base class with files and cuts for IRF analysis

*/
class IRF{
public:
    IRF(std::string summary_root_filename="irf.root");
    ~IRF();

    /// open the input root file(s)
    void open_input_file();

    /// divide a canvas
    void IRF::divideCanvas(TCanvas & c, int nx, int ny, std::string top_title) ;

    const char* angle_cut(int i){
        static char buffer[256];
        sprintf(buffer, "abs(McZDir-%f)<0.1", -0.9+i*0.2);
        return buffer;
    }

    const char* energy_cut(int j){
        static char buffer[256];
        double logecenter=logestart+logedelta*j, ecenter=pow(10, logecenter);
        sprintf(buffer, "abs(McLogEnergy-%f)<%f", logecenter, logedelta/2);
        return buffer;
    }

    double eCenter(int j){
        return pow(10, logestart + j*logedelta);
    }

    const char * hist_name(int i, int j) {
        static char name[8];
        //char * name = new char[8];
        sprintf(name, "h%d%d",j,i);
        return name;
    }

    bool fileExists(){
        TFile f(summary_filename().c_str());
        return f.IsOpen();
    }

    std::string output_file_root(){
        return std::string(::getenv("output_file_root"))+"/";
    }    

    std::string input_file_root(){
        return std::string(::getenv("file_root"))+"/";
    }

    std::string summary_filename(){ return m_summary_filename;}
    std::string input_filename(){return ::getenv("input_file");}

    void set_ymin(double y){ m_ymin=y;}
    void set_ymax(double y){ m_ymax=y;}

    //! optional extra user cut
    void set_user_cut(const TCut & cut){ m_user_cut=cut;}

protected:
    TCut goodCal, goodPSF, goodEvent;

    TCut m_user_cut;

    std::string  m_input_filename;
    std::string  m_summary_filename;

    // define angle bins
    static int angles[5];
    int angle_bins;

    // define energy bins
    double logestart, logedelta;
    int energy_bins;

    // the file
    TFile* m_file;

    // the tree that will be analyzed
    TTree* m_tree;

    // limits used in the draw phase
    double m_ymin, m_ymax;
};
#endif
