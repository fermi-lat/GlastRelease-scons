/** 
      @file makeTot.cxx

  Make sample ToT file in ROOT format using calibRootData classes 
*/
#include <string>
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"   // might not need this one

#include "facilities/Util.h"

#include "commonRootData/idents/TkrId.h"
#include "calibRootData/Tkr/TkrTower.h"
#include "calibRootData/Tkr/Tot.h"

void fillTree(unsigned towerRow, unsigned towerCol, TTree* tree);

int main(int nArg, char** args) {
  std::string outFilename("tot.root");

  if (nArg > 1) outFilename = std::string(args[1]);

  facilities::Util::expandEnvVar(&outFilename);

  TFile* totFile = new TFile(outFilename.c_str(), "RECREATE");
  if (!totFile->IsOpen() ) {
    std::cout << "Could not open file " << outFilename 
              << " for writing " << std::endl;
    return -1;
  }
  else {
    std::cout << "Successfully opened " << outFilename << " for write"
              << std::endl;
  }
  totFile->cd();

  TTree *treeTower00 = new TTree("Tower00", "Tot data for tower 0,0");

  // Do all the serious data fabrication here
  fillTree(0, 0, treeTower00);

  // Add another tree someday

  // finish up
  totFile->cd();
  
  int nBytes = totFile->Write(0, TObject::kOverwrite);

  if (!nBytes) {
    std::cout << "Unable to write file " << outFilename << std::endl; 
    return -1;
  }
  else {
    std::cout << "Wrote file " << outFilename << std::endl;
  }
  totFile->Close();
  delete totFile;

  return 0;
}

void fillTree(unsigned towerRow, unsigned towerCol, TTree* tree) {
  using calibRootData::TkrTower;
  using calibRootData::TotStrip;
  using calibRootData::TotUnilayer;
  using commonRootData::TkrId;

  // Make branch for TkrTower class; 
  TkrTower tower(towerRow, towerCol, TString("hs-tower00"));
  TkrTower* pTower = &tower;
  char* className = "calibRootData::TkrTower";

  TBranch* towerBranch = tree->Branch(className, className, &pTower);

  towerBranch->Fill();  // returns #bytes

  // Make branch for UniLayer class
  //  UInt_t nStrips = 10;  // too lazy to make up 1536 strips
  UInt_t nStrips = 1536;
  TkrId nullId;
  TotUnilayer uni(nullId, nStrips);
  bool sizeOk = uni.consistent();
  if (!sizeOk) {
    std::cout << "Inconsistent sizes.\n  Input size = " << nStrips 
              << std::endl;
    std::cout << "Returned size = " << uni.getNStrips() << std::endl;
  }
  TotUnilayer* pUni = &uni;
  TBranch* uniBranch = tree->Branch("calibRootData::TotUnilayer", 
                                    "calibRootData::TotUnilayer", &pUni);

  // make up a few unilayers.  First is tray 2, top
  uni.setId(TkrId(0,0,2,true,TkrId::eMeasureX));

  TotStrip strip(0, 1.0, 2.0, 3.0, 1.8, 2.0);
  uni.setStrip(strip);

  // strip.init(stripNum, slope, intercept, quad, chi2, df);
  strip.init(1, 1.1, 2.1, 3.1, 1.7, 2.0);
  uni.setStrip(strip);

  strip.init(2, 1.2, 2.2, 3.2, 1.6, 2.0);
  uni.setStrip(strip);

  strip.init(3, 1.3, 2.3, 3.3, 1.5, 2.0);
  uni.setStrip(strip);

  strip.init(4, 1.4, 2.4, 3.4, 1.4, 2.0);
  uni.setStrip(strip);

  strip.init(5, 1.5, 2.5, 3.5, 1.3, 2.0);
  uni.setStrip(strip);

  strip.init(6, 1.6, 2.6, 3.6, 1.2, 2.0);
  uni.setStrip(strip);

  strip.init(7, 1.7, 2.7, 3.7, 1.1, 2.0);
  uni.setStrip(strip);

  strip.init(8, 1.8, 2.8, 3.8, 1.0, 2.0);
  uni.setStrip(strip);

  strip.init(9, 1.8, 2.8, 3.8, 1.0, 2.0);
  uni.setStrip(strip);

  for (unsigned iStrip=10; iStrip < nStrips; iStrip++) {
    strip.init(iStrip, 1.0, 2.0, 3.0, 1.1, 3.0);
    uni.setStrip(strip);
  }
  
  uniBranch->Fill();

  unsigned meas = TkrId::eMeasureY;

  for (unsigned iTray = 3; iTray < 17; iTray++) {
    uni.Clear();

    uni.setId(TkrId(0,0,iTray,false,meas));

    for (unsigned iStrip = 0; iStrip < nStrips; iStrip++) {
      strip.init(iStrip, 2.1, -1.2, 1.2, 1.1, 1.0);
      uni.setStrip(strip);
    }
    uniBranch->Fill();

    uni.Clear();

    uni.setId(TkrId(0,0,iTray,true,meas));

    for (unsigned iStrip = 0; iStrip < nStrips; iStrip++) {
      strip.init(iStrip, 2.1, -1.2, 1.2, 1.1, 1.0);
      uni.setStrip(strip);
    }
    uniBranch->Fill();
    meas = 1 - meas;
  }

}
