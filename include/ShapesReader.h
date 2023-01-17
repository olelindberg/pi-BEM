#ifndef SHAPES_READER_H
#define SHAPES_READER_H

#include "ShapeInput.h"

#include <deal.II/opencascade/boundary_lib.h>
#include <deal.II/opencascade/utilities.h>

#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <Bnd_Box.hxx>

#include <vector>

class ShapesReader
{
public:
  static std::vector<TopoDS_Shape>
  read(const std::filesystem::path &root, const std::vector<ShapeInput> &inputs)
  {
    std::vector<TopoDS_Shape> shapes;
    for (auto &input : inputs)
    {
      auto filename = std::filesystem::path(root).append(input.get_filename());


      std::ifstream file(filename);
      if (file.good())
      {
        std::cout << "Reading CAD shape file: " << filename << std::endl;
        auto shape = dealii::OpenCASCADE::read_IGES(filename, 1e-3);
        shapes.push_back(shape);
      }
      else
      {
        std::cout << filename << std::endl;
        std::cout << "Error" << std::endl;
        // Error handling
      }

      {
        gp_Trsf translate;
        translate.SetTranslation(gp_Vec(input.get_pre_translate_x(),
                                        input.get_pre_translate_y(),
                                        input.get_pre_translate_z()));
        BRepBuilderAPI_Transform translationTransform(translate);
        translationTransform.Perform(shapes.back(), Standard_True);
        shapes.back() = translationTransform.Shape();
      }

      gp_Pnt  origin;
      gp_Trsf scale;
      scale.SetScale(origin, input.get_scale());
      BRepBuilderAPI_Transform scalingTransform(shapes.back(), scale);
      scalingTransform.Perform(shapes.back(), Standard_True);
      shapes.back() = scalingTransform.Shape();

      {
        gp_Trsf rotation;
        rotation.SetRotation(gp_Ax1(origin, gp_Dir(1.0, 0.0, 0.0)), input.get_rotate_x());
        BRepBuilderAPI_Transform transform(shapes.back(), rotation);
        transform.Perform(shapes.back(), Standard_True);
        shapes.back() = transform.Shape();
      }

      {
        gp_Trsf rotation;
        rotation.SetRotation(gp_Ax1(origin, gp_Dir(0.0, 1.0, 0.0)), input.get_rotate_y());
        BRepBuilderAPI_Transform transform(shapes.back(), rotation);
        transform.Perform(shapes.back(), Standard_True);
        shapes.back() = transform.Shape();
      }

      {
        gp_Trsf rotation;
        rotation.SetRotation(gp_Ax1(origin, gp_Dir(0.0, 0.0, 1.0)), input.get_rotate_z());
        BRepBuilderAPI_Transform transform(shapes.back(), rotation);
        transform.Perform(shapes.back(), Standard_True);
        shapes.back() = transform.Shape();
      }

      {
        gp_Trsf translate;
        translate.SetTranslation(gp_Vec(input.get_post_translate_x(),
                                        input.get_post_translate_y(),
                                        input.get_post_translate_z()));
        BRepBuilderAPI_Transform translationTransform(translate);
        translationTransform.Perform(shapes.back(), Standard_True);
        shapes.back() = translationTransform.Shape();
      }
    }
    return shapes;
  }
};

#endif // SHAPES_READER_H
