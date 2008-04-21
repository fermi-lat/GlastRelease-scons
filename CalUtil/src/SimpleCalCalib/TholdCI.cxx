/** @file
    @author Zachary Fewtrell
*/

// LOCAL INCLUDES
#include "CalUtil/SimpleCalCalib/TholdCI.h"

// GLAST INCLUDES

// EXTLIB INCLUDES

// STD INCLUDES
#include <sstream>
#include <string>
#include <fstream>

using namespace CalUtil;
using namespace std;

namespace CalUtil {
  const short TholdCI::INVALID_THOLD = -5000;

  TholdCI::TholdCI() :
    m_fleThresh(FaceIdx::N_VALS, INVALID_THOLD),
    m_fheThresh(FaceIdx::N_VALS, INVALID_THOLD),
    m_lacThresh(FaceIdx::N_VALS, INVALID_THOLD),
    m_uldThresh(RngIdx::N_VALS, INVALID_THOLD),
    m_ped(RngIdx::N_VALS, INVALID_THOLD)
  {
  }

  void TholdCI::readTXT(const string &filename) {
    // open file
    ifstream infile(filename.c_str());
    
    if (!infile.is_open())
      throw runtime_error(string("Unable to open " + filename));
    
    // loop through each line in file
    string line;
    unsigned short twr,lyr,col,face;
    float lac, fle, fhe;
    CalVec<RngNum, float> uld, ped;
    while (infile.good()) {
      getline(infile, line);
      if (infile.fail()) break; // bad get

      // check for comments
      if (line[0] == ';')
        continue;

      istringstream istrm(line);

      istrm >> twr >> lyr >> col >> face 
            >> lac >> fle >> fhe 
            >> uld[0] >> uld[1] >> uld[2] >> uld[3]
            >> ped[0] >> ped[1] >> ped[2] >> ped[3];

      const FaceIdx faceIdx(twr,
                            LyrNum(lyr),
                            col,
                            // the parenthesis are _very_ important.
                            // see Scott Meyers "Effective STL" Item 6
                            (FaceNum(face)));
      
      m_fleThresh[faceIdx] = fle;
      m_fheThresh[faceIdx] = fhe;
      m_lacThresh[faceIdx] = lac;

      for (RngNum rng; rng.isValid(); rng++) {
        const RngIdx rngIdx(faceIdx,rng);
        m_uldThresh[rngIdx] = uld[rng];
        m_ped[rngIdx] = ped[rng];
      } // rng loop
    } // line loop
  } // readTXT()
} // namespace CalUtil
  
