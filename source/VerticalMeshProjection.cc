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

// #  include <BRepAdaptor_CompCurve.hxx>
// #  include <BRepAdaptor_Curve.hxx>
// #  include <BRepAdaptor_HCompCurve.hxx>
// #  include <BRepAdaptor_HCurve.hxx>
// #  include <BRepTools.hxx>
// #  include <BRep_Tool.hxx>
// #  include <GCPnts_AbscissaPoint.hxx>
// #  include <ShapeAnalysis_Curve.hxx>
// #  include <ShapeAnalysis_Surface.hxx>
// #  include <Standard_Version.hxx>
// #  include <TopoDS.hxx>

#include <deal.II/opencascade/utilities.h>

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_HCompCurve.hxx>
#include <BRepAdaptor_HCurve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgo_Section.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <Extrema_ExtAlgo.hxx>
#include <Extrema_ExtFlag.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_Line.hxx>
#include <IntCurvesFace_ShapeIntersector.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>

#include "../include/VerticalMeshProjection.h"
#include "../include/my_utilities.h"
#include <limits>

DEAL_II_NAMESPACE_OPEN

using namespace OpenCASCADE;


template <int dim, int spacedim>
VerticalMeshProjection<dim, spacedim>::VerticalMeshProjection(TopoDS_Shape &sh,
                                                              double        reference_level,
                                                              const double  tolerance)
  : sh(sh)
  , _reference_level(reference_level)
  , tolerance(tolerance)
{
  TopExp_Explorer exp;
  for (exp.Init(sh, TopAbs_FACE); exp.More(); exp.Next())
  {
    TopoDS_Face face = TopoDS::Face(exp.Current());

    Bnd_Box bnd_box;
    BRepBndLib::Add(face, bnd_box);

    bvh_vec_t          point_min(bnd_box.CornerMin().X(), bnd_box.CornerMin().Y());
    bvh_vec_t          point_max(bnd_box.CornerMax().X(), bnd_box.CornerMax().Y());
    BVH_Box<double, 2> bvh_box(point_min, point_max);

    _bvh_boxset.Add(face, bvh_box);
  }
  _bvh_boxset.Build();
}

template <int dim, int spacedim>
std::unique_ptr<Manifold<dim, spacedim>> VerticalMeshProjection<dim, spacedim>::clone() const
{
  return std::unique_ptr<Manifold<dim, spacedim>>(
    new VerticalMeshProjection<dim, spacedim>(sh, _reference_level, tolerance));
}

void assign_to_bvh_surface_selector(const BVH_SurfaceSelector &bvh_surface_selector,
                                    const bvh_boxset_t &       bvh_boxset,
                                    const bvh_vec_t &          point)
{
  BVH_SurfaceSelector &bvhss = const_cast<BVH_SurfaceSelector &>(bvh_surface_selector);
  bvhss.SetBVHSet(const_cast<bvh_boxset_t *>(&bvh_boxset));
  bvhss.set_point(point);
}

template <int dim, int spacedim>
Point<spacedim>
VerticalMeshProjection<dim, spacedim>::project_to_manifold(const ArrayView<const Point<spacedim>> &,
                                                           const Point<spacedim> &candidate) const
{
  dealii::Point<3> point(candidate[0], candidate[1], _reference_level);

  if (std::fabs(candidate[0] - 1100.0) < 1.0e-1)
    std::cout << "upper bnd\n";

  BVH_SurfaceSelector bvh_surface_selector;
  assign_to_bvh_surface_selector(bvh_surface_selector,
                                 _bvh_boxset,
                                 bvh_vec_t(candidate[0], candidate[1]));
  if (bvh_surface_selector.Select() > 0)
  {
    Handle(Geom_Line) line2 =
      GC_MakeLine(gp_Pnt(candidate[0], candidate[1], 1000.0), gp_Dir(0.0, 0.0, -1.0));
    for (auto &surface : bvh_surface_selector.get_surfaces())
    {
      GeomAPI_IntCS geom_intersect;
      geom_intersect.Perform(line2, surface);
      if (geom_intersect.IsDone() && geom_intersect.NbPoints() > 0)
      {
        point[2] = geom_intersect.Point(1).Z();
      }
    }
  }

  if (point[2] > 10.0)
    point[2] = 10.0;
  return point;
}

template class VerticalMeshProjection<2, 3>;

DEAL_II_NAMESPACE_CLOSE
