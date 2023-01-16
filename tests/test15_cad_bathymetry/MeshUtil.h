#ifndef MESH_UTIL_H
#define MESH_UTIL_H

#include <deal.II/base/geometry_info.h>

#include <vector>

class MeshUtil
{
public:
  static void
  project(dealii::Triangulation<2, 3> &mesh, unsigned int manifold_id)
  {
    for (auto cell = mesh.begin_active(); cell != mesh.end(); ++cell)
    {
      if (cell->manifold_id() == manifold_id)
      {
        std::vector<dealii::Point<3>> vertices(dealii::GeometryInfo<2>::vertices_per_cell);
        for (unsigned int i = 0; i < dealii::GeometryInfo<2>::vertices_per_cell; ++i)
          vertices[i] = cell->vertex(i);

        for (unsigned int v = 0; v < dealii::GeometryInfo<2>::vertices_per_cell; ++v)
        {
          cell->vertex(v) = cell->get_manifold().project_to_manifold(vertices, cell->vertex(v));
        }
      }
    }
  }
};

#endif // MESH_UTIL_H
