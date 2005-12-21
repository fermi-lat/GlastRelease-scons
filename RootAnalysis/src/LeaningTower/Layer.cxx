#include "Layer.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

Layer::Layer(TString name, float pz, float py, float px,
             float rz, float ry, float rx) {
    SetName(name);
    Z = pz;
    Y = py;
    X = px;
    rotZ = rz;
    rotY = ry;
    rotX = rx;
    dX = dY = dZ = dRotX = dRotY = dRotZ = 0.0;

    for ( int i=0; i<4; i++ ) {
        float xcoord1 = GetCoordinate(i*384);
        float xcoord2 = GetCoordinate((i+1)*384-1);
        LadderLine[i] = new TLine(xcoord1, Z, xcoord2, Z);
    }
    LayerLine  = new TLine(GetCoordinate(0), Z, GetCoordinate(1535), Z);
    LayerLabel = new TText(400, Z, fName);
    LayerLabel->SetTextAlign(02);
    LayerLabel->SetTextSize(0.02f);
}

//////////////////////////////////////////////////
Layer::~Layer() {
    //    std::cout<<"~Layer()"<<std::endl;
    delete LayerLine;
    delete LayerLabel;
    for ( int i=0; i<4; i++ )
        delete LadderLine[i];
}

//////////////////////////////////////////////////
float Layer::GetCoordinate(int StripNumber) {
    bool DEBUG = false;
    // Check that the strip number is within the allowed range.
    if ( StripNumber < 0 ||  StripNumber > 1535 ) {
        std::cout << "WARNING: the strip number must be in the range "
                  << "[0:1535]." << std::endl;
        std::cout << "GetCoordinate() returning -1.0" << std::endl;
        return -1.0;
    }

    /*
    float pos = INACTIVEBORDERWIDTH + STRIPPITCH * (StripNumber+0.5)
        + ( LADDERGAP + 2 * INACTIVEBORDERWIDTH )
        * static_cast<float>(StripNumber/384);
        */
    float pos = GetLadderXmin(StripNumber/384) + INACTIVEBORDERWIDTH
        + STRIPPITCH * ( StripNumber%384 + 0.5 );

    if ( DEBUG )
        std::cout << "GetCoordinate() returning " << pos <<
            "(strip number = " << StripNumber << ")" << std::endl;
    return pos;
}

/*
bool Layer::checkActiveArea(float ParallelCoordinate, float NormalCoordinate,
                            float BorderWidth) {
  bool DEBUG = false;
  // Define some constants (to be moved out of the function?)
  // All dimension in mm.
  // Loop over the ladders.

  // THIS HERE IS NOT PRECISE, BUT I DON'T WANT TO MESS WITH IT!
  for (int i=0; i<4; i++){
    if ((NormalCoordinate > (INACTIVEBORDERWIDTH + (SIWAFERSIDE + LADDERGAP)*i
                             + BorderWidth)) &&
	(NormalCoordinate < (SIWAFERSIDE*(i+1) + LADDERGAP*i
                             - BorderWidth))){
      for (int j=0; j<4; j++){
	if ((ParallelCoordinate > (INACTIVEBORDERWIDTH + SIWAFERSIDE*j
        + BorderWidth)) &&
	    (ParallelCoordinate < (SIWAFERSIDE*(j+1) - BorderWidth))){
	  if (DEBUG){
	    std::cout << "checkActiveArea(): ladder " << i << ", wafer " << j
                      << " hit." << std::endl;
	  }
	  return 1;
	}
      }
    }
  }
  if (DEBUG){
    std::cout << "checkActiveArea(): no hit on any of the wafers." << std::endl;
  }
  return 0;
}
*/

float Layer::activeAreaDist(const float x, const float y) const {
    // returns the distance of the hit to an active area.
    // if negative, the hit lies in an active area.
    // Translations of the planes are accounted for in Recon::TkrAlignmentSvc.
    // ToDo:
    //   Also correct for rotations!

    float dx = 1000000.f;
    for ( int i=0; i<4; ++i ) {
        // x dist to ladder center
        const float dist = std::abs(GetLadderCenX(i)-x);
        if ( dist < dx )
            dx = dist;
    }
    float dy = 1000000.f;
    for ( int i=0; i<4; ++i ) {
        // y dist to wafer center
        const float dist = std::abs(GetWaferCenY(i)-y);
        if ( dist < dy )
            dy = dist;
    }

    dx -= SIWAFERACTIVESIDE2;
    dy -= SIWAFERACTIVESIDE2;
    //      |   
    //      |  if both distances are positiv, the distance is the sqrt(dx^2+dy^2)
    //      |  otherwise it's the maximum of both
    //------+
    //
    //
    return ( dx>0 && dy>0 ) ? sqrt(dx*dx+dy*dy) : std::max(dx, dy);
}

void Layer::DrawLayer() {
    for ( int i=0; i<4; i++ ) {
        //      LadderLine[i]->SetLineStyle(0);
        LadderLine[i]->Draw();
    }
    LayerLabel->Draw();
}

void Layer::DrawGhostLayer()
{
  LayerLine->SetLineColor(17);
  LayerLine->Draw();
}

void Layer::SetTree(TFile *file) {
    TString LayerTreeName = "Layer";
    LayerTreeName += GetName();
    LayerTree = (TTree*) file->Get(LayerTreeName);
    LayerTree->SetBranchAddress("TkrNumHits",&TkrNumHits);
    LayerTree->SetBranchAddress("TkrHits",TkrHits);
    LayerTree->SetBranchAddress("TriggerReq0",&TriggerReq0);
    LayerTree->SetBranchAddress("TriggerReq1",&TriggerReq1);
    LayerTree->SetBranchAddress("ToT0",&ToT0);
    LayerTree->SetBranchAddress("ToT1",&ToT1);
}

std::string Layer::SaveGeometry() const {
    std::stringstream s;
    s.setf(std::ios_base::fixed);
    s.precision(3);
    s << fName << ' ' << Z+dZ << ' ' << Y+dY << ' ' << X+dX;
    s.precision(2);
    s << ' ' << rotZ+dRotZ << ' ' << rotY+dRotY << ' ' << rotX+dRotX;
    return s.str();
}

ClassImp(Layer)
