//  class galacticHist
//  Author:  Theodore Hierath
//  This class contains a histogram and graphing class for
//  energies.
//  The bins are numbered from 0 to number of bins - 1;

#include "galacticHist.h"
#include <iostream>
#include <iomanip>
#include <fstream>

galacticHist::galacticHist(){}

galacticHist::~galacticHist(void)
{}
//////////////////////TEST CODE////////////////////////////
void galacticHist::test(){

    std::ofstream out_file("graph.cxx", std::ios::app);

    out_file.clear();
    //out_file << 
           // "{\n"
      //      "   gROOT->Reset();\n";

      out_file << 
"{\n"
//"#include <iostream.h>\n"
//"#include <stdio.h>\n"
//"#include <fcntl.h>\n"

"  FILE *fp;\n"
"  float ptmp,p[20];\n"
"  int i, iline=0;\n";
out_file << 
"  ntuple = new TNtuple(" << '"' << "ntuple" << '"' << "," << '"' << "NTUPLE" << '"' << "," << '"' <<"x:y:z" << '"' << ");\n";

out_file <<
"  fp = fopen(" << '"' << "./data.dat" << '"' << "," << '"' << "r" << '"' << ");\n"

"  while ( fscanf(fp," << '"' << "%f" << '"' << ",&ptmp) != EOF ){\n"
"    p[i++]=ptmp;\n"
"    if (i==3){\n"
"      i=0; \n"
"      iline++;\n"
"      ntuple->Fill(p[0],p[1],p[2]); \n"
"    }\n"
"  }\n"


"  ntuple.Draw(" << '"' << "z:y:x" << '"' 
//<< "," << '"' << '"'
//<< "," << '"' << "HIST" << '"' 
<< ");\n"

//"  ntuple.Draw(" << '"' << "x>>h1" << '"' << "," << '"' << '"' << "," << '"' << "h1" << '"' << ");\n"
//"  h1.Draw("
//<< '"' << "CONT" << '"'
//<<");\n"

"}\n";

out_file.close();

system("root -l graph.cxx");


}