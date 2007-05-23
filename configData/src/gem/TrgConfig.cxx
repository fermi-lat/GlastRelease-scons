//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      A top level GEM configuration class 
//
// Environment:
//      Software developed for GLAST.
//
// Author List:
//      Martin Kocian
//
// Copyright Information:
//      Copyright (C) 2005      Stanford Linear Accelerator Center
//
//---------------------------------------------------------------------------
#include "configData/gem/TrgConfig.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>


TrgConfig::TrgConfig(){}

void TrgConfig::printContrigurator(std::ostream& out) const{
  int width[]={7,9,10,9,4,7,7,4,4,10,7,9,7,8};
  
  out<<setiosflags(std::ios::left);
  out<<std::setfill(' ')<<std::setw(width[0])<<"Engine";
  out<<std::setfill(' ')<<std::setw(width[1])<<"External";
  out<<std::setfill(' ')<<std::setw(width[2])<<"Solicited";
  out<<std::setfill(' ')<<std::setw(width[3])<<"Periodic";
  out<<std::setfill(' ')<<std::setw(width[4])<<"CNO";
  out<<std::setfill(' ')<<std::setw(width[5])<<"CAL-HI";
  out<<std::setfill(' ')<<std::setw(width[6])<<"CAL-LO";
  out<<std::setfill(' ')<<std::setw(width[7])<<"TKR";
  out<<std::setfill(' ')<<std::setw(width[8])<<"ROI";
  out<<std::setfill(' ')<<std::setw(width[9])<<"Zero-supp";
  out<<std::setfill(' ')<<std::setw(width[10])<<"Ranges";
  out<<std::setfill(' ')<<std::setw(width[11])<<"Prescale";
  out<<std::setfill(' ')<<std::setw(width[12])<<"Marker";
  out<<std::setfill(' ')<<std::setw(width[13])<<"Inhibit";
  out<<std::endl;
  for (int i=0;i<16;i++){
    out<<std::setfill(' ')<<std::setw(width[0])<<i;
    std::vector<unsigned long> condvec=lut()->conditionsInEngine(i);
    for (int cond=7;cond>=0;cond--){
      bool zero=false;
      bool one=false;
      for (unsigned j=0;j<condvec.size();j++){
	if(condvec[j]&(1<<cond))one=true;
	else zero=true;
      }
      char en=' ';
      if (zero&&one)en='x';
      else if (zero)en='0';
      else if (one)en='1';
      out<<std::setfill(' ')<<std::setw(width[8-cond])<<en;
    }
    std::string zs;
    if (trgEngine()->zeroSuppression(i))zs="enabled";
    else zs="disabled";
    out<<std::setfill(' ')<<std::setw(width[9])<<zs;
    if (trgEngine()->fourRangeReadout(i))zs="4-rgn";
    else zs="1-rng";
    out<<std::setfill(' ')<<std::setw(width[10])<<zs;
    int prescale=trgEngine()->prescale(i);
    out<<std::setfill(' ')<<std::setw(width[11])<<prescale;
    int marker=trgEngine()->marker(i);
    out<<std::setfill(' ')<<std::setw(width[11])<<marker;
    if (trgEngine()->inhibited(i))zs="inhibit";
    else zs="no";
    out<<std::setfill(' ')<<std::setw(width[10])<<zs;
    out<<std::endl;
  }
}


// Reset the cached and output values
void TrgConfig::reset() {
  _twp.reset();
  _tpt.reset();
  _lut.reset();
  _tev.reset();
  _tdv.reset();
  _configuration.reset();
  _roi.reset();  
}
std::ostream& operator <<(std::ostream &os, const TrgConfig& tl){
  os<<"GEM configuration: ";
  os<<*tl.configuration()<<std::endl;
  os<<*tl.windowParams()<<std::endl;
  os<<*tl.periodicTrigger()<<std::endl;
  os<<*tl.lut()<<std::endl;
  os<<*tl.trgEngine()<<std::endl;
  os<<*tl.disabledChannels()<<std::endl;
  os<<*tl.roi()<<std::endl;
  return os;
}

// Attach this value to a TTree
void TrgConfig::makeBranches(TTree& tree, const std::string& prefix) const {
  _twp.makeBranch(tree,prefix);
  _tpt.makeBranch(tree,prefix);
  _lut.makeBranch(tree,prefix);
  _tev.makeBranch(tree,prefix);
  _tdv.makeBranch(tree,prefix);
  _configuration.makeBranch(tree,prefix);
  _roi.makeBranch(tree,prefix);  
 
}

void TrgConfig::attach(TTree& tree, const std::string& prefix) const {
  _twp.attach(tree,prefix);
  _tpt.attach(tree,prefix);
  _lut.attach(tree,prefix);
  _tev.attach(tree,prefix);
  _tdv.attach(tree,prefix);
  _configuration.attach(tree,prefix);
  _roi.attach(tree,prefix);  
}
