#ifndef HEIGHT_RATIO_REFINEMENT_H
#define HEIGHT_RATIO_REFINEMENT_H

#include "GridRefinement.h"

#include <iostream>
#include <vector>

class HeightRatioRefinement : public GridRefinement
{
public:
  HeightRatioRefinement(dealii::ConditionalOStream pcout, int itermax_, double height_ratio_max_)
    : GridRefinement(pcout)
    , itermax(itermax_)
    , height_ratio_max(height_ratio_max_)
  {}

  virtual ~HeightRatioRefinement()
  {}

  virtual void refine(dealii::Triangulation<2, 3> &tria) override;

private:
  int    itermax          = 0;
  double height_ratio_max = 0.0;
};

#endif // HEIGHT_RATIO_REFINEMENT_H