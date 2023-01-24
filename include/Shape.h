#ifndef SHAPE_H
#define SHAPE_H

#include <BRepBuilderAPI_Transform.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Ax1.hxx>

class Shape
{
public:
  // Set position in world coordinates
  void set_position(double valx, double valy, double valz)
  {
    double translate_x = valx - _position_x;
    double translate_y = valy - _position_y;
    double angle_z     = _rotate_z;

    double translate_x_local = std::cos(angle_z) * translate_x - std::sin(angle_z) * translate_y;
    double translate_y_local = std::sin(angle_z) * translate_x + std::cos(angle_z) * translate_y;

    {
      gp_Trsf translate;
      translate.SetTranslation(gp_Vec(translate_x_local, translate_y_local, 0.0));
      BRepBuilderAPI_Transform transform(translate);
      transform.Perform(shape);
      shape = transform.Shape();
    }

    _position_x = valx;
    _position_y = valy;
    _position_z = valz;
  }

  void set_rotation(double valx, double valy, double valz)
  {
    gp_Pnt origin;
    {
      gp_Trsf rotation;
      rotation.SetRotation(gp_Ax1(origin, gp_Dir(0.0, 0.0, 1.0)), valz - _rotate_z);
      BRepBuilderAPI_Transform transform(rotation);
      transform.Perform(shape, Standard_True);
      shape = transform.Shape();
    }

    _rotate_z = valz;
  }

  TopoDS_Shape shape;

private:
  double _position_x = 0.0;
  double _position_y = 0.0;
  double _position_z = 0.0;

  double _rotate_x = 0.0;
  double _rotate_y = 0.0;
  double _rotate_z = 0.0;
};

#endif // SHAPE_H