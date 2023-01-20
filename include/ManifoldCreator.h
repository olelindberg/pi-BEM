#ifndef MANIFOLD_CREATOR_H
#define MANIFOLD_CREATOR_H

#include "IPhysicalDomain.h"
#include "Shape.h"
#include "VerticalMeshProjection.h"
#include "my_manifold_lib.h"

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/opencascade/boundary_lib.h>
#include <deal.II/opencascade/utilities.h>

#include <memory>

class ManifoldCreator
{
public:
  static std::vector<std::shared_ptr<dealii::Manifold<2, 3>>>
  make(const std::vector<ShapeInput> &inputs,
       std::vector<Shape> &           shapes,
       dealii::Triangulation<2, 3> &  mesh)
  {
    std::vector<std::shared_ptr<dealii::Manifold<2, 3>>> manifolds;

    auto shape = shapes.begin();
    for (auto &input : inputs)
    {
      if (input.get_mesh_projection() == "directional")
      {
        std::cout << "creating directional projection" << std::endl;
        auto manifold = std::make_shared<dealii::VerticalMeshProjection<2, 3>>(
          shape->shape,
          input.get_pre_translate_z(),
          1.0 * dealii::OpenCASCADE::get_shape_tolerance(shape->shape));
        manifolds.push_back(manifold);
      }
      else if (input.get_mesh_projection() == "normal")
      {
        std::cout << "creating normal projection" << std::endl;
        auto manifold = std::make_shared<dealii::MyNormalToMeshProjectionManifold<2, 3>>(
          shape->shape, 1.0 * dealii::OpenCASCADE::get_shape_tolerance(shape->shape));
        manifolds.push_back(manifold);
      }

      else if (input.get_mesh_projection() == "arch_length")
      {
        std::cout << "creating arch length projection" << std::endl;
        auto manifold =
          std::make_shared<dealii::OpenCASCADE::ArclengthProjectionLineManifold<2, 3>>(
            shape->shape, 1.0 * dealii::OpenCASCADE::get_shape_tolerance(shape->shape));
        manifolds.push_back(manifold);
      }

      mesh.set_manifold(input.get_mesh_element_id(), *manifolds.back());

      ++shape;
    }
    return manifolds;
  }
};
#endif // MANIFOLD_CREATOR_H