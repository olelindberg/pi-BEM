#ifndef SHAPES_READER_H
#define SHAPES_READER_H

#include "Shape.h"
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
  static std::vector<Shape> read(const std::filesystem::path &  root,
                                 const std::vector<ShapeInput> &inputs)
  {
    std::vector<Shape> shapes;
    for (auto &input : inputs)
    {
      auto filename = std::filesystem::path(root).append(input.get_filename());

      Shape shape;

      std::ifstream file(filename);
      if (file.good())
      {
        std::cout << "Reading CAD shape file: " << filename << std::endl;
        shape.shape = dealii::OpenCASCADE::read_IGES(filename, 1e-3);
      }
      else
      {
        std::cout << filename << std::endl;
        std::cout << "Error" << std::endl;
        // Error handling
      }

      shape.set_position(input.get_pre_translate_x(),
                         input.get_pre_translate_y(),
                         input.get_pre_translate_z());

      shape.set_rotation(input.get_rotate_x(), input.get_rotate_y(), input.get_rotate_z());

      shapes.push_back(shape);
    }
    return shapes;
  }
};

#endif // SHAPES_READER_H
