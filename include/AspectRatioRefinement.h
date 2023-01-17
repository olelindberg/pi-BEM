#ifndef ASPECT_RATIO_REFINEMENT_H
#define ASPECT_RATIO_REFINEMENT_H

#include "GridRefinement.h"

#include <iostream>
#include <vector>

class AspectRatioRefinement : public GridRefinement
{
public:
  AspectRatioRefinement(dealii::ConditionalOStream pcout,
                        const std::vector<int> &   manifold_id_,
                        int                        itermax_,
                        double                     aspect_ratio_max_)
    : GridRefinement(pcout)
    , manifold_id(manifold_id_)
    , itermax(itermax_)
    , aspect_ratio_max(aspect_ratio_max_)
  {}

  virtual ~AspectRatioRefinement()
  {}

  virtual void refine(dealii::Triangulation<2, 3> &tria) override;

private:
  std::vector<int> manifold_id;
  int              itermax          = 0;
  double           aspect_ratio_max = 0.0;
};

#endif // ASPECT_RATIO_REFINEMENT_H