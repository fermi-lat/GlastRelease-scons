// $Header$

// Include files
// Gaudi system includes
#include "idents/CalXtalId.h"

#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"
#include "CalUtil/CalArray.h"

#include <vector>

using namespace CalUtil;
using namespace idents;

bool testCalVec(){
  cout << "testCalVec" << endl;

  // creats vector of size 8, which is
  // enought to test functionality
  CalVec<XtalRng, unsigned> testVec;

  testVec.resize(XtalRng::N_VALS);
  
  if (testVec.size() != (unsigned)XtalRng::N_VALS)
    return false;

  for (XtalRng idx; idx.isValid(); idx++) {
    testVec[idx] = idx.val();
  }

  for (XtalRng idx; idx.isValid(); idx++) {
    if (testVec[idx] != idx.val())
      return false;
    
    if (testVec.find(idx.val()) - testVec.begin() 
        != (int)idx.val())
      return false;
  }

  testVec.fill(3);
  for (XtalRng idx; idx.isValid(); idx++) {
    if (testVec[idx] != 3)
      return false;
  }

  testVec.clear();
  if (testVec.size() != (unsigned)0) 
    return false;

  return true;
}

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
    
    if ((testArr.find(idx.val()) - testArr.begin()) 
        != (int)idx.val())
      return false;
  }

  testArr.fill(3);
  for (XtalRng idx; idx.isValid(); idx++) {
    if (testArr[idx] != 3)
      return false;
  }

  return true;
}

bool testCalDefs(){

  cout << "testCalDefs" << endl;

  //-- TEST CALDEFS CONVERSION ROUTINES --//
  // try to test all conversions & constructors for
  // consistency.
  for (short tRow=0; tRow < 4 ; tRow++)
    for (short tCol=0; tCol < 4; tCol++) {

      TwrNum twr(tRow, tCol);
      if (tRow != twr.getRow()) return false;
      if (tCol != twr.getCol()) return false;
      
      for (DirNum dir; dir.isValid(); dir++)
        for (GCRCNum gcrc; gcrc.isValid(); gcrc++) {

          LyrNum lyr(dir, gcrc);
          if (dir != lyr.getDir()) return false;

          for (ColNum col; col.isValid(); col++) {
            CalXtalId xtalId(twr,lyr,col);

            XtalIdx xtalIdx(twr,lyr,col);
            if (xtalId != xtalIdx.getCalXtalId()) return false;
            if (twr != xtalIdx.getTwr()) return false;
            if (lyr != xtalIdx.getLyr()) return false;
            if (col != xtalIdx.getCol()) return false;

            // alt-constructor(s)
            if (XtalIdx(xtalId) != xtalIdx) return false;

            for (FaceNum face; face.isValid(); face++) {
              FaceIdx faceIdx(twr,lyr,col,face);
              CalXtalId faceId(twr,lyr,col,(CalXtalId::XtalFace)face);
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

                XtalDiode xDiode(face,diode);
                if (diode != xDiode.getDiode()) return false;
                if (face  != xDiode.getFace()) return false;

                DiodeIdx diodeIdx(twr,lyr,col,face,diode);
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

                  RngNum rng(diode,thx);
                  if (diode != rng.getDiode()) return false;

                  XtalRng xRng(face, rng);
                  if (face != xRng.getFace()) return false;
                  if (rng != xRng.getRng()) return false;
                  if (xDiode != xRng.getXtalDiode()) return false;

                  RngIdx rngIdx(twr,lyr,col,face,rng);
                  CalXtalId rngId(twr, lyr, col, face.val(), rng.val());
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

int main()
{
  bool status = true;
  
  if (!testCalDefs())  status = false;
  if (!testCalVec())   status = false;
  if (!testCalArray()) status = false;
  
  if (!status) {
    cout << "You are a failure!" << endl;
    return -1;
  }
  
  cout << "You are a success!" << endl;
  return 0;
}


