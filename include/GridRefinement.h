#ifndef GRID_REFINEMENT_H
#define GRID_REFINEMENT_H


#include "IGridRefinement.h"

#include <deal.II/base/conditional_ostream.h>


class GridRefinement : public IGridRefinement
{
public:
  GridRefinement(dealii::ConditionalOStream pcout, bool verbose = false)
    : _pcout(pcout)
    , _verbose(verbose)
  {}
  virtual ~GridRefinement()
  {}

  virtual void refine(dealii::Triangulation<2, 3> &){};

protected:
  dealii::ConditionalOStream _pcout;
  bool                       _verbose = false;
};

#endif // GRID_REFINEMENT_H
