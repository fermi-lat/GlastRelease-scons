#include "detModel/Utilities/GDDmatrix.h"

void GDDtranslate(GDDvector* t, GDDvector* v)
{
  v->x = v->x + t->x;
  v->y = v->y + t->y;
  v->z = v->z + t->z;
};

/// 0 = x, 1 = y, 2 = z
void GDDrotate(int dir, double theta, GDDvector* v)
{
  GDDvector* temp = new GDDvector(v->x, v->y, v->z);

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
