#ifndef GRID_REFINEMENT_H
#define GRID_REFINEMENT_H


#include "IGridRefinement.h"

#include <deal.II/base/conditional_ostream.h>


class GridRefinement : public IGridRefinement
{
public:
  GridRefinement(dealii::ConditionalOStream pcout)
    : _pcout(pcout)
  {}
  virtual ~GridRefinement()
  {}

  virtual void
  refine(dealii::Triangulation<2, 3> &tria, const std::vector<TopoDS_Shape> &cad_surfaces){};

protected:
  dealii::ConditionalOStream _pcout;

private:
};

#endif // GRID_REFINEMENT_H
