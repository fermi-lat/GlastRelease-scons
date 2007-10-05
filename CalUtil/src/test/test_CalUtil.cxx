// $Header$

/** @file 
    @author Zach Fewtrell
    
    CalUtil test app. 
    Simple non-gaudi main() method runs test method for each CalUtil module
    
    @return non-zero to os on failure
 */

// LOCAL INLUDES
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"
#include "CalUtil/CalArray.h"

// GLAST INCLUDES
#include "idents/CalXtalId.h"

// EXTLIB INCLUDES

// STD INCLUDES
#include <vector>

using namespace CalUtil;
using namespace idents;

/// test CalVec class
bool testCalVec(){
  cout << "testCalVec" << endl;

  // creats vector of size 8, which is
  // enought to test functionality
  CalVec<XtalRng, unsigned> testVec;

  testVec.resize(XtalRng::N_VALS);
  
  if (testVec.size() != (unsigned)XtalRng::N_VALS)
    return false;

  /// set each value to value of index
  for (XtalRng idx; idx.isValid(); idx++)
    testVec[idx] = idx.val();


  for (XtalRng idx; idx.isValid(); idx++) {
    if (testVec[idx] != idx.val())
      return false;
  }


  testVec.clear();
  if (testVec.size() != (unsigned)0) 
    return false;

  return true;
}

/// test CalArray class
bool testCalArray(){

  cout << "testCalArray" << endl;

  // creats array of size 8, which is
  // enough to test functionality
  CalArray<XtalRng, unsigned> testArr;

  if (testArr.size() != (unsigned)XtalRng::N_VALS)
    return false;

  for (XtalRng idx; idx.isValid(); idx++) {
    testArr[idx] = idx.val();
  }

  for (XtalRng idx; idx.isValid(); idx++) {
    if (testArr[idx] != idx.val())
      return false;
    
  }

  return true;
}

/// \brief test suite of CalDefs classes
///
/// Basically loop through all Cal components & check constructors & converters over all types and all crystals
bool testCalDefs(){

  cout << "testCalDefs" << endl;

  //-- TEST CALDEFS CONVERSION ROUTINES --//
  // try to test all conversions & constructors for
  // consistency.
  for (short tRow=0; tRow < 4 ; tRow++)
    for (short tCol=0; tCol < 4; tCol++) {

      const TwrNum twr(tRow, tCol);
      if (tRow != twr.getRow()) return false;
      if (tCol != twr.getCol()) return false;
      
      for (DirNum dir; dir.isValid(); dir++)
        for (GCRCNum gcrc; gcrc.isValid(); gcrc++) {

          const LyrNum lyr(dir, gcrc);
          if (dir != lyr.getDir()) return false;

          for (ColNum col; col.isValid(); col++) {
            const CalXtalId xtalId(twr.val(),
                                   lyr.val(),
                                   col.val());

            const XtalIdx xtalIdx(twr,lyr,col);
            if (xtalId != xtalIdx.getCalXtalId()) return false;
            if (twr != xtalIdx.getTwr()) return false;
            if (lyr != xtalIdx.getLyr()) return false;
            if (col != xtalIdx.getCol()) return false;

            // alt-constructor(s)
            if (XtalIdx(xtalId) != xtalIdx) return false;

            for (FaceNum face; face.isValid(); face++) {
              const FaceIdx faceIdx(twr,lyr,col,face);
              const CalXtalId faceId(twr.val(),
                                     lyr.val(),
                                     col.val(),
                                     (CalXtalId::XtalFace)face.val());

              if (faceId != faceIdx.getCalXtalId()) return false;
              if (twr != faceIdx.getTwr()) return false;
              if (lyr != faceIdx.getLyr()) return false;
              if (col != faceIdx.getCol()) return false;
              if (face.val() != faceIdx.getFace().val()) return false;
              if (xtalIdx != faceIdx.getXtalIdx()) return false;

              // alt-constructor(s)
              if (FaceIdx(faceId) != faceIdx) return false;
              if (FaceIdx(xtalIdx, face) != faceIdx) return false;

              
              for (DiodeNum diode; diode.isValid(); diode++) {
                if (RngNum(diode,THX8).val() != diode.getX8Rng().val()) return false;
                if (RngNum(diode,THX1) != diode.getX1Rng()) return false;

                const XtalDiode xDiode(face,diode);
                if (diode != xDiode.getDiode()) return false;
                if (face  != xDiode.getFace()) return false;

                const DiodeIdx diodeIdx(twr,lyr,col,face,diode);
                if (twr != diodeIdx.getTwr()) return false;
                if (lyr != diodeIdx.getLyr()) return false;
                if (col != diodeIdx.getCol()) return false;
                if (face != diodeIdx.getFace()) return false;
                if (diode != diodeIdx.getDiode()) return false;
                if (xtalIdx != diodeIdx.getXtalIdx()) return false;
                if (faceIdx != diodeIdx.getFaceIdx()) return false;

                //alt-constructor(s)
                if (DiodeIdx(xtalIdx, face, diode) != diodeIdx) return false;
                if (DiodeIdx(xtalIdx, xDiode) != diodeIdx) return false;
                if (DiodeIdx(faceIdx, diode) != diodeIdx) return false;
                
              
                
                for (THXNum thx; thx.isValid(); thx++) {

                  const RngNum rng(diode,thx);
                  if (diode != rng.getDiode()) return false;

                  const XtalRng xRng(face, rng);
                  if (face != xRng.getFace()) return false;
                  if (rng != xRng.getRng()) return false;
                  if (xDiode != xRng.getXtalDiode()) return false;

                  const RngIdx rngIdx(twr,lyr,col,face,rng);
                  const CalXtalId rngId(twr.val(), 
                                        lyr.val(), 
                                        col.val(), 
                                        face.val(), 
                                        rng.val());
                  if (twr != rngIdx.getTwr()) return false;
                  if (lyr != rngIdx.getLyr()) return false;
                  if (col != rngIdx.getCol()) return false;
                  if (face != rngIdx.getFace()) return false;
                  if (rng != rngIdx.getRng()) return false;
                  if (xtalIdx != rngIdx.getXtalIdx()) return false;
                  if (faceIdx != rngIdx.getFaceIdx()) return false;
                  
                  //alt-constructor(s)
                  if (RngIdx(xtalIdx, face, rng) != rngIdx) return false;
                  if (RngIdx(xtalIdx, xRng) != rngIdx) return false;
                  if (RngIdx(faceIdx, rng) != rngIdx) return false;


                } // thx/rng
              } // diode
            } // face
          } // col
        } // gcrc
    } // twr

  return true;
}

/// main method tests each CalUtil component in sequence
/// \return non-zero to OS on failure
int main()
{
  bool status = true;
  
  if (!testCalDefs())  status = false;
  if (!testCalVec())   status = false;
  if (!testCalArray()) status = false;
  
  if (!status) {
    cout << "CalUtil unit test failed!" << endl;
    return -1;
  }
  
  cout << "You are a success!" << endl;
  return 0;
}


