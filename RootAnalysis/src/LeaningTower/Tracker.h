#include "TMap.h"
#include "TString.h"

class Tracker
{
  
 public:
  Tracker();
  void loadGeometry(TString filename="src/LeaningTower/geometry/stack2geometry.txt");
  void Display();
  inline  TMap* GetGeometry() {return myGeometry;}
  inline void IsTower(bool tower) {TOWER=tower;}
 private:
  bool TOWER;
  TMap* myGeometry;
};

