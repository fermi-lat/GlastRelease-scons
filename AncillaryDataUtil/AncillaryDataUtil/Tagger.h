#ifndef TAGGER_HH
#define TAGGER_HH

// ROOT classes.
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TVectorT.h>
#include <TF1.h>
#include <iostream>
#include <fstream>

// Our own classes.
#include "AncillaryDataEvent/TaggerModule.h"
#include "AncillaryDataUtil/AncillaryDataServer.h"

class Tagger

{
 public:
 
  Tagger();
  ~Tagger();

  TaggerModule  *getModule(int moduleId);
  TaggerLayer   *getLayer(int moduleId, int layerId);
  TaggerChannel *getChannel(int moduleId, int layerId, int channelId);
  void readPedestals(std::string pedestalFilePath);
  void readEvent(AncillaryDataEvent *event);
  void processData(AncillaryDataEvent *event);
  void resetData();
  void writeOutputRootFile(std::string outputRootFilePath);
 
 private:
  
  bool m_DEBUG;
  TaggerModule* m_modules[N_MODULES];
  std::vector<double> m_trackX;
  std::vector<double> m_trackY;
  std::vector<double> m_trackZ;
  TGraph   *m_trackXZ;
  TGraph   *m_trackYZ;
};
#endif
