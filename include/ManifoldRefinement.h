#ifndef MANIFOLD_REFINEMENT_H
#define MANIFOLD_REFINEMENT_H

#include "GridRefinement.h"

class ManifoldRefinement : public GridRefinement
{
public:
  ManifoldRefinement(dealii::ConditionalOStream pcout, unsigned int manifoldId_, int levels_)
    : GridRefinement(pcout)
    , manifoldId(manifoldId_)
    , levels(levels_)
  {}
  virtual ~ManifoldRefinement()
  {}

  virtual void
  refine(dealii::Triangulation<2, 3> &tria, const std::vector<TopoDS_Shape> &cad_surfaces) override;

private:
  unsigned int manifoldId = 0;
  int          levels     = 0;
};

#endif // MANIFOLD_REFINEMENT_H
