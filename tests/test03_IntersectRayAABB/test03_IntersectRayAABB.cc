#include "my_utilities.h"
#include <iostream>

int
main ()
{
  dealii::Point<3>     origin (-0.642229, -0.642229, -0.288675);
  dealii::Tensor<1, 3> direct;
  direct[0] = 0.672922;
  direct[1] = 0.674851;
  direct[2] = 0.302905;

  for (int testid = 0; testid < 2; ++testid)
  {
    dealii::Point<3> *corner_min;
    dealii::Point<3> *corner_max;
    switch (testid)
    {
    case 0:
      corner_min = new dealii::Point<3> (-0.707403, -1.00029, -0.707369);
      corner_max = new dealii::Point<3> (0.707403, -0.577054, 0.000296207);
      break;
    case 1:
      corner_min = new dealii::Point<3> (-0.707403, 0.577054, -0.707369);
      corner_max = new dealii::Point<3> (0.707403, 1.00029, 0.000296436);
      break;
    default:
      break;
    }

    double           tmin;
    dealii::Point<3> point;
    std::cout << "test " << testid;
    if (IntersectRayAABB (origin, direct, *corner_min, *corner_max, tmin, point))
      std::cout << ": Intersection\n";
    else
      std::cout << ": No intersection\n";
  }
  return 0;
}
