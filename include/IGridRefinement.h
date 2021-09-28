#ifndef IGRID_REFINEMENT_H
#define IGRID_REFINEMENT_H

#include <deal.II/grid/tria.h>

#include <TopoDS_Shape.hxx>


class IGridRefinement
{
public:
  IGridRefinement()
  {}
  virtual ~IGridRefinement()
  {}

  virtual void
  refine(dealii::Triangulation<2, 3> &tria, const std::vector<TopoDS_Shape> &cad_surfaces) = 0;
};

#endif // IGRID_REFINEMENT_H
