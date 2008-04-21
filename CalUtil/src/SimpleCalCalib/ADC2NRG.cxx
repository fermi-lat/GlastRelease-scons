// $Header$

/** @file
    @author Zachary Fewtrell
*/

// LOCAL INCLUDES
#include "CalUtil/SimpleCalCalib/ADC2NRG.h"

// GLAST INCLUDES

// EXTLIB INCLUDES

// STD INCLUDES
#include <fstream>
#include <sstream>

using namespace CalUtil;
using namespace std;


namespace CalUtil {

  const float ADC2NRG::INVALID_ADC2NRG = -5000;

  ADC2NRG::ADC2NRG() :
    m_adc2nrg(RngIdx::N_VALS, INVALID_ADC2NRG)
  {
  }

  void ADC2NRG::writeTXT(const string &filename) const {
    ofstream outfile(filename.c_str());

    // output header info as comment
    outfile << ";twr lyr col face rng adc2nrg" << endl;

    if (!outfile.is_open())
      throw runtime_error(string("Unable to open " + filename));

    for (RngIdx rngIdx; rngIdx.isValid(); rngIdx++)
      outfile << rngIdx.getTwr().val()
              << " " << rngIdx.getLyr().val()
              << " " << rngIdx.getCol().val()
              << " " << rngIdx.getFace().val()
              << " " << rngIdx.getRng().val()
              << " " << m_adc2nrg[rngIdx]
              << endl;
  }

  void ADC2NRG::readTXT(const string &filename) {
    ifstream infile(filename.c_str());

    if (!infile.is_open())
      throw runtime_error(string("Unable to open " + filename));

    string line;
    while (infile.good()) {
      float adc2nrg;
      unsigned short twr;
      unsigned short lyr;
      unsigned short col;
      unsigned short face;
      unsigned short rng;

      getline(infile, line);
      if (infile.fail()) break; // bad get

      // check for comments
      if (line[0] == ';')
        continue;

      istringstream istrm(line);

      istrm >> twr
            >> lyr
            >> col
            >> face
            >> rng
            >> adc2nrg;

      const RngIdx rngIdx(twr,
                          LyrNum(lyr),
                          col,
                          FaceNum((idents::CalXtalId::XtalFace)face),
                          rng);

      m_adc2nrg[rngIdx] = adc2nrg;
    }
  }

}; // namespace CalUtil
