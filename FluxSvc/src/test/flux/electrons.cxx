#include <iostream>
#include <stdlib.h>
#include <utility>

// CLHEP
#include <CLHEP/Random/JamesRandom.h>


#include "../CrElectron.h"

typedef  double G4double;

const int BIN_NUM = 1024;
const int LOOP    = int(1e+6);

const G4double energy_min = 0.1;
const G4double energy_max = 100.0;


int main(int argc, char** argv)
{
	using namespace std;
	const double GeV=1, rad=1;
  int hist[BIN_NUM];
  for (int i = 0; i < BIN_NUM; i++) hist[i] = 0;

  CrElectron src;

  HepRandomEngine* engine = new HepJamesRandom;

  int loop = (argc > 1) ? (int)atof(*(argv + argc - 1)) : LOOP;
  {for (int i = 0; i < loop; i++){
    src.selectComponent(engine);
    G4double energy = src.energySrc(engine) / GeV;
    std::pair<double,double> dir = src.dir(energy, engine);
      // Notice: dir consists of cos(zenith_angle) and azimuth [rad]
    G4double theta  = acos(dir.first);  // [0, pi]
    G4double phi    = dir.second * rad;
    hist[int(BIN_NUM * (energy - energy_min) / (energy_max - energy_min))]++;
    if (i % 10000 == 1) cerr << '\015' << i << ": " << energy << "...";
  }}
  cerr << "\n";

  {for (int i = 0; i < BIN_NUM; i++){
    cout << ((i + 0.5)/BIN_NUM) * (energy_max - energy_min) + energy_min << " "
         << hist[i] << "\n";
  }}

  delete engine;

  return 0;
}
