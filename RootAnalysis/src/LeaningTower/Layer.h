#include "TLine.h"
#include "TText.h"

#include "TTree.h"
#include "TFile.h"

class Layer : public TObject
{
 public:
  Layer();
  ~Layer();
  inline void SetLayerName(TString name){Name = name;}
  inline  TString GetLayerName() {return Name;}
  void SetHeight(double height);
  inline double GetHeight(){return Z;}

  double GetCoordinate(int stripId);
  bool checkActiveArea(double ParallelCoordinate, double NormalCoordinate, double BorderWidth);
  void DrawLayer();
  void DrawGhostLayer();
  inline bool IsX(){return Name.Contains("X");}
  inline bool IsY(){return Name.Contains("Y");}

  ClassDef(Layer,1);

  TString Name;
  double ShiftX;
  double ShiftY;
  double Z;
  TLine *LayerLine;
  TLine *LadderLine[4];
  TTree *LayerTree;
  TText *LayerLabel;
  void SetTree(TFile *file);
  void GetEvent(int event);

  int TkrNumHits;
  int TkrHits[128];

  inline void AddHitInActiveArea()
    {
      HitsInActiveArea++;
    }

  inline void AddMissedHit()
    {
      MissedHits++;
    }

  inline double GetInefficiency()
    {
      if(HitsInActiveArea)
	return (double) MissedHits/HitsInActiveArea;
    }
  inline double GetEfficiency()
    {
      return 1.0-GetInefficiency();
    }

  Bool_t GetTriggerReq(Bool_t side) { return side ? TriggerReq1 : TriggerReq0; }
  
 private:
  double EDGE_WIDTH,WAFER_WIDTH,STRIP_PITCH,LADDER_SEPARATION;
  int HitsInActiveArea,MissedHits;
  Bool_t TriggerReq0;
  Bool_t TriggerReq1;

};
