#ifndef GDDBOUNDINGBOX_H
#define GDDBOUNDINGBOX_H

#include <iostream>
#define MAX(a,b) (a>b)?a:b
#define MIN(a,b) (a<b)?a:b

/**
 * This class represent the bounding box (axis aligned) of a volume 
 * @author R.Giannitrapani
 */
class GDDboundingBox {
 public:
  enum axisDir{xDir, yDir, zDir};

 public:
  /// A constructors 
  GDDboundingBox():xDim(0),yDim(0),zDim(0){;};
  /// A constructors with dimension
  GDDboundingBox(double x, double y, double z):xDim(x),yDim(y),zDim(z){;};
  /**
   * This methods sets and gets the bounding box dimensions 
   */
  void setXDim(double x){xDim = x;};
  void setYDim(double y){yDim = y;};
  void setZDim(double z){zDim = z;};
  double getXDim(){return xDim;};
  double getYDim(){return yDim;};
  double getZDim(){return zDim;};

  void setDirDim(axisDir dir, double a);
  double getDirDim(axisDir dir);

  /** 
   * This method update the dimensions for a rotation of an angle
   * around one of the axis
   */
  void rotate(double xAngle, double yAngle, double zAngle);
  void rotateDir(axisDir dir, double a);

  bool exists();

  /** 
   * This methods update the dimensions for the merging of another
   * bounding box translated by (x,y,z) wrt the center of the
   * original bounding box
   */
  void merge(GDDboundingBox* b, double x, double y, double z);

 private:
  /// The x dimension of the bounding box
  double xDim;
  /// The y dimension of the bounding box
  double yDim;
  /// The z dimension of the bounding box
  double zDim;
};
#endif //GDDBOUNDINGBOX_H

