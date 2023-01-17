#ifndef BOX_REFINEMENT_H
#define BOX_REFINEMENT_H

#include "GridRefinement.h"
#include <Bnd_Box.hxx>

class BoxRefinement : public GridRefinement
{
public:
  BoxRefinement(dealii::ConditionalOStream pcout,
                const Bnd_Box &            box,
                double                     aspectRatioMax,
                double                     cellSizeMin,
                unsigned int               manifold_id_,
                int                        levels_);
  virtual ~BoxRefinement(){};

  virtual void refine(dealii::Triangulation<2, 3> &tria) override;

private:
  Bnd_Box _box;
  double  _aspectRatioMax = 0.0;
  double  _cellSizeMin    = 0.0;

  unsigned int manifold_id = 0;
  int          levels      = 0;
};
#endif // BOX_REFINEMENT_H
