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

DEAL_II_NAMESPACE_OPEN

using namespace OpenCASCADE;


template <int dim, int spacedim>
VerticalMeshProjection<dim, spacedim>::VerticalMeshProjection(TopoDS_Shape &sh,
                                                              double        reference_level,
                                                              const double  tolerance)
  : sh(sh)
  , _reference_level(reference_level)
  , tolerance(tolerance)
{}

template <int dim, int spacedim>
std::unique_ptr<Manifold<dim, spacedim>> VerticalMeshProjection<dim, spacedim>::clone() const
{
  return std::unique_ptr<Manifold<dim, spacedim>>(
    new VerticalMeshProjection<dim, spacedim>(sh, _reference_level, tolerance));
}


template <int dim, int spacedim>
Point<spacedim>
VerticalMeshProjection<dim, spacedim>::project_to_manifold(const ArrayView<const Point<spacedim>> &,
                                                           const Point<spacedim> &candidate) const
{
  Tensor<1, 3> average_normal;
  average_normal[0] = 0.0;
  average_normal[1] = 0.0;
  average_normal[2] = -1.0;

  dealii::Point<3> point;
  if (!my_line_intersection(sh, candidate, average_normal, tolerance, point))
  {
    point = dealii::Point<3>(candidate[0], candidate[1], _reference_level);
  }
  return point;
}

template class VerticalMeshProjection<2, 3>;

DEAL_II_NAMESPACE_CLOSE
