// $Header$

/** @file
    @author Zachary Fewtrell
*/

// LOCAL INCLUDES
#include "CalUtil/SimpleCalCalib/CalDAC.h"

// GLAST INCLUDES

// EXTLIB INCLUDES

// STD INCLUDES
#include <fstream>
#include <sstream>

using namespace CalUtil;
using namespace std;


namespace CalUtil {

  const int CalDAC::INVALID_DAC = -5000;

  CalDAC::CalDAC() :
    m_DACs(FaceIdx::N_VALS, INVALID_DAC)
  {
  }

  void CalDAC::writeTXT(const string &filename) const {
    ofstream outfile(filename.c_str());

    // output header info as comment
    outfile << ";twr lyr col face DAC" << endl;

    if (!outfile.is_open())
      throw runtime_error(string("Unable to open " + filename));

    for (FaceIdx faceIdx; faceIdx.isValid(); faceIdx++)
      outfile << faceIdx.getTwr().val()
              << " " << faceIdx.getLyr().val()
              << " " << faceIdx.getCol().val()
              << " " << faceIdx.getFace().val()
              << " " << m_DACs[faceIdx]
              << endl;
  }

  void CalDAC::readTXT(const string &filename) {
    ifstream infile(filename.c_str());

    if (!infile.is_open())
      throw runtime_error(string("Unable to open " + filename));

    string line;
    while (infile.good()) {
      unsigned short DAC;
      unsigned short twr;
      unsigned short lyr;
      unsigned short col;
      unsigned short face;

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
            >> DAC;

      const FaceIdx faceIdx(twr,
                            LyrNum(lyr),
                            col,
                            FaceNum((idents::CalXtalId::XtalFace)face));

      m_DACs[faceIdx]   = DAC;
    }
  }

}; // namespace CalUtil
