#include "detModel/Utilities/Matrix.h"
#include <math.h>

namespace detModel{

  void Translate(Vector* t, Vector* v)
  {
    v->x = v->x + t->x;
    v->y = v->y + t->y;
    v->z = v->z + t->z;
  };

  void Equate(Vector* t, Vector* v)
  {
    t->x = v->x;
    t->y = v->y;
    t->z = v->z;
  };

  /// 0 = x, 1 = y, 2 = z
  void Rotate(int dir, double theta, Vector* v)
  {
    Vector* temp = new Vector(v->x, v->y, v->z);

    theta = theta*GDDPI/180;

    switch(dir){
    case 0:
      {
	v->x = temp->x;
	v->y = temp->y*cos(theta) - temp->z*sin(theta);
	v->z = temp->y*sin(theta) + temp->z*cos(theta);
      }
      break;

    case 1:
      {
	v->x = temp->x*cos(theta) + temp->z*sin(theta);
	v->y = temp->y;
	v->z = temp->z*cos(theta) - temp->x*sin(theta);
      }
      break;

    case 2:
      {
	v->x = temp->x*cos(theta) - temp->y*sin(theta);
	v->y = temp->x*sin(theta) + temp->y*cos(theta);
	v->z = temp->z;
      }
      break;
    }

    delete temp;

  }
}
