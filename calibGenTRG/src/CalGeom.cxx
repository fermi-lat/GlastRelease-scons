/** @file CalGeom.cxx

@author Zachary Fewtrell

\brief function definitions for CalGeom.h

$Header$
*/

// LOCAL INCLUDES
#include "CalGeom.h"

// GLAST INCLUDES

// EXTLIB INCLUDES

// STD INCLUDES

using namespace CalUtil;

    CalUtil::XtalIdx pos2Xtal(const Vec3D &pos) {
      // tower 'columns' numbered in X dimension
      short twrCol = (short)(pos.x()/twrPitch + 2);


      if (twrCol < 0 || twrCol > 3) return INVALID_XTAL;

      // tower 'rows' numbered in Y dimension
      short twrRow = (short)(pos.y()/twrPitch + 2);
      if (twrRow < 0 || twrRow > 3) return INVALID_XTAL;
      TwrNum twr(twrRow,
                 twrCol);

      float zDiff  = CalTopZ - pos.z();
      LyrNum lyr( (short)(zDiff / cellVertPitch));

      // check that we are in a valid xtal
      if (!lyr.isValid())
        return INVALID_XTAL;

      // check that we are not in region between this xtal and next
      float  zOffset =  zDiff-lyr.val()*cellVertPitch;
      if ( zOffset > CsIHeight)
        return INVALID_XTAL;

      ColNum col;
      if (lyr.getDir() == X_DIR) {
        float twrCtrY = (twrRow-1.5)*twrPitch;
        col = (short)((pos.y() - twrCtrY)/cellHorPitch + 6);

        // check that we are not in xtal gap
        float xtalCtr = cellHorPitch*((float)col.val()-5.5)+twrCtrY;
        if (abs(pos.y() - xtalCtr) > CsIWidth/2)
          return INVALID_XTAL;
      } else {
        // Y_DIR
        float twrCtrX = (twrCol-1.5)*twrPitch;
        col = (short)((pos.x() - twrCtrX)/cellHorPitch + 6);
        // check that we are not in xtal gap
        float xtalCtr = cellHorPitch*((float)col.val()-5.5)+twrCtrX;
        if (abs(pos.x() - xtalCtr) > CsIWidth/2)
          return INVALID_XTAL;
      }

      if (!col.isValid())
        return INVALID_XTAL;

      return XtalIdx(twr, lyr, col);
    }

    Vec3D xtalCtrPos(const CalUtil::XtalIdx xtalIdx) {
      const LyrNum lyr     = xtalIdx.getLyr();
      const ColNum col     = xtalIdx.getCol();
      const TwrNum twr     = xtalIdx.getTwr();

      const float  z       = lyr0ZCtr - cellVertPitch*lyr.val();
      float  x, y;
      const float  twrCtrY = ((float)twr.getRow()-1.5)*twrPitch;
      const float  twrCtrX = ((float)twr.getCol()-1.5)*twrPitch;


      if (lyr.getDir() == X_DIR) {
        x = twrCtrX;
        y = twrCtrY + cellHorPitch*((float)col.val()-5.5);
      } else {
        y = twrCtrY;
        x = twrCtrX + cellHorPitch*((float)col.val()-5.5);
      }

      return Vec3D(x, y, z);
    }
