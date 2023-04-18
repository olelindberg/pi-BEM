#ifndef AREA_REFINEMENT_H
#define AREA_REFINEMENT_H

#include "GridRefinement.h"

#include <Bnd_Box.hxx>

class AreaRefinement : public GridRefinement
{
public:
  AreaRefinement(dealii::ConditionalOStream pcout,
                 const std::vector<int>    &manifold_id_,
                 double                     aspectRatioMax,
                 int                        levels_,
                 double                     areaMax,
                 bool                       verbose);
  virtual ~AreaRefinement(){};

  virtual void refine(dealii::Triangulation<2, 3> &tria) override;

private:
  std::vector<int> manifold_id;

  int levels = 0;

  double _aspectRatioMax = 0.0;
  double _areaMax        = 0.0;
};
#endif // AREA_REFINEMENT_H
