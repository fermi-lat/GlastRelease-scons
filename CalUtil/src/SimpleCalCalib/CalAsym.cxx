// $Header$

/** @file
    @author Zachary Fewtrell
*/

// LOCAL INCLUDES
#include "CalUtil/SimpleCalCalib/CalAsym.h"

// GLAST INCLUDES

// EXTLIB INCLUDES
#include "TH2S.h"

// STD INCLUDES
#include <fstream>
#include <sstream>

using namespace std;
using namespace CalUtil;
  
namespace CalUtil {
  CalAsym::CalAsym(const unsigned short nSlicesPerXtal) :
    nSlicesPerXtal(nSlicesPerXtal)
  {
  }

  void CalAsym::writeTXT(const std::string &filename) const {
    ofstream outfile(filename.c_str());

    if (!outfile.is_open())
      throw runtime_error(string("Unable to open " + filename));

    // output header info as comment
    outfile << ";twr lyr col pos_diode neg_diode asym error" << endl;

    // PER XTAL LOOP
    for (XtalIdx xtalIdx; xtalIdx.isValid(); xtalIdx++) {
      const TwrNum twr (xtalIdx.getTwr());
      const LyrNum lyr (xtalIdx.getLyr());
      const ColNum col (xtalIdx.getCol());
      for (AsymType asymType; asymType.isValid(); asymType++)
        // per point along curve
        for (unsigned short i = 0; i < nSlicesPerXtal; i++) {
          // skip empty channels
          if (m_asym[asymType][xtalIdx].size() != nSlicesPerXtal)
            continue;

          outfile << twr.val()
                  << " " << lyr.val()
                  << " " << col.val()
                  << " " << asymType.getDiode(POS_FACE).val()
                  << " " << asymType.getDiode(NEG_FACE).val()
                  << " " << m_asym[asymType][xtalIdx][i]
                  << " " << m_asymErr[asymType][xtalIdx][i]
                  << endl;
        }
    }
  }

  void CalAsym::readTXT(const std::string &filename) {
    unsigned short twr, lyr, col, pdiode, ndiode;
    float asym, sig;

    // open file
    ifstream infile(filename.c_str());

    if (!infile.is_open())
      throw runtime_error(string("Unable to open " + filename));

    // loop through each line in file
    string line;
    while (infile.good()) {
      getline(infile, line);
      if (infile.fail()) break; // bad get

      // check for comments
      if (line[0] == ';')
        continue;

      istringstream istrm(line);

      // get lyr, col (xtalId)
      istrm >> twr >> lyr >> col >> pdiode >> ndiode >> asym >> sig;

      const XtalIdx xtalIdx(twr,
                            LyrNum(lyr),
                            col);
      
      // the parenthesis are _very_ important.
      // see Scott Meyers "Effective STL" Item 6
      AsymType asymType((DiodeNum(pdiode)),(DiodeNum(ndiode)));
      
      m_asym[asymType][xtalIdx].push_back(asym);
      m_asymErr[asymType][xtalIdx].push_back(sig);
    }
  }

  void CalAsym::genSplines() {
    for (DiodeNum diode; diode.isValid(); diode++) {
      fill(m_a2pSplines[diode].begin(), m_a2pSplines[diode].end(), static_cast<TSpline3*>(0));
      fill(m_p2aSplines[diode].begin(), m_p2aSplines[diode].end(), static_cast<TSpline3*>(0));
    }

    // create position (Y-axis) array
    // linearly extrapolate for 1st and last points (+2 points)
    double *pos = new double[nSlicesPerXtal+2];
    for (unsigned short i = 0; i < nSlicesPerXtal+2; i++)
      pos[i] = i + 0.5; // (center of the column)
    double *asym = new double[nSlicesPerXtal+2];

    // PER XTAL LOOP
    for (DiodeNum diode; diode.isValid(); diode++) {
      const AsymType asymType(diode,
                              diode);

      for (XtalIdx xtalIdx; xtalIdx.isValid(); xtalIdx++) {
        // skip empty channels
        if (m_asym[asymType][xtalIdx].size() != nSlicesPerXtal)
          continue;

        // copy asym vector into middle of array
        vector<float> &asymVec(m_asym[asymType][xtalIdx]);
        copy(asymVec.begin(), asymVec.end(), asym+1);

        // extrapolate 1st & last points
        asym[0] = 2*asym[1] - asym[2];
        asym[nSlicesPerXtal+1]     = 2*asym[nSlicesPerXtal]-asym[nSlicesPerXtal-1];

        {
          //generate splinename
          ostringstream name;
          name << "asym2pos_" << xtalIdx.val() << "_" << diode.val();

          // create spline object
          TSpline3     *mySpline = new TSpline3(name.str().c_str(),
                                                asym, pos, nSlicesPerXtal+2);
          mySpline->SetName(name.str().c_str());

          m_a2pSplines[diode][xtalIdx] = mySpline;
        }

        // create inverse spline object
        {
          //generate splinename
          ostringstream name;
          name << "pos2asym_" << xtalIdx.val() << "_" << diode.val();

          // create spline object
          TSpline3     *mySpline = new TSpline3(name.str().c_str(),
                                                pos, asym, nSlicesPerXtal+2);
          mySpline->SetName(name.str().c_str());

          m_p2aSplines[diode][xtalIdx] = mySpline;
        }
      }
    }

        delete [] pos;
        delete [] asym;
  }

}; // namespace CalUtil
