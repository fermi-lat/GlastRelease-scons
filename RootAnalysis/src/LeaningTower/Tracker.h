#include "TMap.h"
#include "TString.h"

class Tracker
{
  
 public:
  Tracker();
  void loadGeometry(TString filename="src/LeaningTower/geometry/stack2geometry.txt");
  void Display();
  inline  TMap* GetGeometry() {return myGeometry;}
 private:
  TMap* myGeometry;
};

