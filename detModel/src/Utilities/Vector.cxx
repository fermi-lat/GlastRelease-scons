#include "detModel/Utilities/Vector.h"
#include <math.h>

namespace detModel{
  
  void Vector::Add(Vector* t)
  {
    x += t->x;
    y += t->y;
    z += t->z;
  };

  void Vector::Equate(Vector* v)
  {
    x = v->x;
    y = v->y;
    z = v->z;
  };

  /// 0 = x, 1 = y, 2 = z
  void Vector::Rotate(int dir, double theta)
  {
    Vector* temp = new Vector(x, y, z);

    theta = theta*GDDPI/180;

    switch(dir){
    case 0:
      {
	x = temp->x;
	y = temp->y*cos(theta) - temp->z*sin(theta);
	z = temp->y*sin(theta) + temp->z*cos(theta);
      }
      break;

    case 1:
      {
	x = temp->x*cos(theta) + temp->z*sin(theta);
	y = temp->y;
	z = temp->z*cos(theta) - temp->x*sin(theta);
      }
      break;

    case 2:
      {
	x = temp->x*cos(theta) - temp->y*sin(theta);
	y = temp->x*sin(theta) + temp->y*cos(theta);
	z = temp->z;
      }
      break;
    }

    delete temp;

  }


  void Vector::Multiply(float a)
  {
    x = x*a;
    y = y*a;
    z = z*a;
  }

}
