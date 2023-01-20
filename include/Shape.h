#ifndef SHAPE_H
#define SHAPE_H

#include <BRepBuilderAPI_Transform.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Ax1.hxx>

class Shape
{
public:
  void set_position(double valx, double valy, double valz)
  {
    {
      gp_Pnt  origin;
      gp_Trsf rotation;
      rotation.SetRotation(gp_Ax1(origin, gp_Dir(0.0, 0.0, 1.0)), -_rotate_z);
      BRepBuilderAPI_Transform transform(rotation);
      transform.Perform(shape, Standard_True);
      shape = transform.Shape();
    }
    {
      gp_Trsf translate;
      translate.SetTranslation(gp_Vec(valx - _position_x, valy - _position_y, valz - _position_z));
      BRepBuilderAPI_Transform transform(translate);
      transform.Perform(shape);
      shape = transform.Shape();
    }
    {
      gp_Pnt  origin;
      gp_Trsf rotation;
      rotation.SetRotation(gp_Ax1(origin, gp_Dir(0.0, 0.0, 1.0)), _rotate_z);
      BRepBuilderAPI_Transform transform(rotation);
      transform.Perform(shape, Standard_True);
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

    {
      gp_Trsf rotation;
      rotation.SetRotation(gp_Ax1(origin, gp_Dir(0.0, 1.0, 0.0)), valy - _rotate_y);
      BRepBuilderAPI_Transform transform(rotation);
      transform.Perform(shape, Standard_True);
      shape = transform.Shape();
    }

    {
      gp_Trsf rotation;
      rotation.SetRotation(gp_Ax1(origin, gp_Dir(1.0, 0.0, 0.0)), valx - _rotate_x);
      BRepBuilderAPI_Transform transform(rotation);
      transform.Perform(shape, Standard_True);
      shape = transform.Shape();
    }

    _rotate_x = valx;
    _rotate_y = valy;
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