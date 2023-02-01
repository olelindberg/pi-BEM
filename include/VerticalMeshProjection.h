// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#ifndef VERTICAL_MESH_PROJECTION_H
#define VERTICAL_MESH_PROJECTION_H

#include "../include/BVH_SurfaceSelector.h"

#include <deal.II/base/config.h>

#include <deal.II/grid/manifold.h>
#include <deal.II/opencascade/utilities.h>
#include <unordered_map>


DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
class VerticalMeshProjection : public FlatManifold<dim, spacedim>
{
public:
  VerticalMeshProjection(TopoDS_Shape &sh, double reference_level, const double tolerance = 1e-7);
  virtual std::unique_ptr<Manifold<dim, spacedim>> clone() const override;
  virtual Point<spacedim>
  project_to_manifold(const ArrayView<const Point<spacedim>> &surrounding_points,
                      const Point<spacedim> &                 candidate) const override;

protected:
  mutable std::unordered_map<std::size_t, Tensor<1, spacedim>> projections_cache;
  TopoDS_Shape &                                               sh;
  double                                                       _reference_level = 0.0;
  const double                                                 tolerance;
  bvh_boxset_t                                                 _bvh_boxset;


private:
};


/*@}*/

DEAL_II_NAMESPACE_CLOSE


#endif // dealii_occ_manifold_lib_h
