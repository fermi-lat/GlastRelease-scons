#include "Layer.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>



//////////////////////////////////////////////////
Layer::Layer()
{
  ShiftX=0;
  ShiftY=0;
 
  EDGE_WIDTH        = 0.1000;
  STRIP_PITCH       = 0.0228;
  LADDER_SEPARATION = 0.0200; 
  WAFER_WIDTH       = 8.9500;

  MissedHits=0;
  HitsInActiveArea=0;
}
//////////////////////////////////////////////////
Layer::~Layer()
{
  std::cout<<"~Layer()"<<std::endl;
  delete LayerLine;
  delete LayerLabel;
  for(int i=0;i<4;i++)
    {
      delete LadderLine[i];
    }
}
//////////////////////////////////////////////////
void Layer::SetHeight(double height)
{
  Z=height;
  for(int i=0;i<4;i++)
    {
      double xcoord1 = GetCoordinate(i*384);
      double xcoord2 = GetCoordinate((i+1)*384-1);
      
      
      LadderLine[i] = new TLine(xcoord1, height, xcoord2, height);
    }
  LayerLine = new TLine(GetCoordinate(0), height, GetCoordinate(1535),height);
  
}
//////////////////////////////////////////////////

double Layer::GetCoordinate(int StripNumber)
{
  bool DEBUG = false;
  // Check that the strip number is within the allowed range.
  if ((StripNumber < 0) || (StripNumber > 1535)){
    std::cout << "WARNING: the strip number must be included in the range [0:1535]." << std::endl;
    std::cout << "GetCoordinate() returning -1.0" << std::endl;
    return -1.0;
  }
  // Define some constants (to be moved out of the function?)
  // All dimension in cm.
  
  double Coordinate = EDGE_WIDTH + STRIP_PITCH*StripNumber +
    (LADDER_SEPARATION + 2*EDGE_WIDTH - STRIP_PITCH)*(int)(StripNumber/384);
  if (DEBUG)
    {
      std::cout << "GetCoordinate() returning " << Coordinate <<
	"(strip number = " << StripNumber << ")" << std::endl;
    }
  return Coordinate;
}

bool Layer::checkActiveArea(double ParallelCoordinate, double NormalCoordinate, double BorderWidth)
{
  bool DEBUG = false;
  // Define some constants (to be moved out of the function?)
  // All dimension in cm.
  // Loop over the ladders.
  for (int i=0; i<4; i++){
    if ((NormalCoordinate > (EDGE_WIDTH + (WAFER_WIDTH + LADDER_SEPARATION)*i + BorderWidth)) &&
	(NormalCoordinate < (WAFER_WIDTH*(i+1) + LADDER_SEPARATION*i - BorderWidth))){
      for (int j=0; j<4; j++){
	if ((ParallelCoordinate > (EDGE_WIDTH + WAFER_WIDTH*j + BorderWidth)) &&
	    (ParallelCoordinate < (WAFER_WIDTH*(j+1) - BorderWidth))){
	  if (DEBUG){
	    std::cout << "checkActiveArea(): ladder " << i << ", wafer " << j << " hit." << std::endl;
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

void Layer::DrawLayer()
{
  for(int i=0;i<4;i++)
    {
      //      LadderLine[i]->SetLineStyle(0);
      LadderLine[i]->Draw();
    }
  LayerLabel = new TText(38, Z, Name);
  LayerLabel->SetTextAlign(02);
  LayerLabel->SetTextSize(0.02);
  LayerLabel->Draw();
  
}

void Layer::DrawGhostLayer()
{
  LayerLine->SetLineColor(17);
  LayerLine->Draw();
}
void Layer::SetTree(TFile *file)
{
  TString LayerTreeName = "Layer";
  LayerTreeName+=Name;
  LayerTree = (TTree*) file->Get(LayerTreeName);
  LayerTree->SetBranchAddress("TkrNumHits",&TkrNumHits);
  LayerTree->SetBranchAddress("TkrHits",TkrHits);
  LayerTree->SetBranchAddress("TriggerReq0",&TriggerReq0);
  LayerTree->SetBranchAddress("TriggerReq1",&TriggerReq1);
}

void Layer::GetEvent(int event)
{
  LayerTree->GetEntry(event);
}

