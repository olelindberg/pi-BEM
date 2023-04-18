#ifndef DISTANCE_REFINEMENT_H
#define DISTANCE_REFINEMENT_H

#include "GridRefinement.h"

#include <Bnd_Box.hxx>

class DistanceRatioRefinement : public GridRefinement
{
public:
  DistanceRatioRefinement(dealii::ConditionalOStream pcout,
                          const Bnd_Box             &box,
                          double                     aspectRatioMax,
                          double                     area_min,
                          unsigned int               manifold_id_,
                          int                        levels_,
                          double                     distanceRatioMax,
                          bool                       verbose = false);
  virtual ~DistanceRatioRefinement(){};

  virtual void refine(dealii::Triangulation<2, 3> &tria) override;

private:
  Bnd_Box _box;
  double  _boxLengthMax = 0.0;

  double _aspectRatioMax = 0.0;
  double _area_min       = 0.0;

  unsigned int manifold_id       = 0;
  int          levels            = 0;
  double       _distanceRatioMax = 0.0;
};
#endif // DISTANCE_REFINEMENT_H
