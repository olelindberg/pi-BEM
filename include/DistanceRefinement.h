#ifndef DISTANCE_REFINEMENT_H
#define DISTANCE_REFINEMENT_H

#include "GridRefinement.h"

#include <Bnd_Box.hxx>

class DistanceRefinement : public GridRefinement
{
public:
  DistanceRefinement(dealii::ConditionalOStream pcout,
                     const Bnd_Box &            box,
                     unsigned int               manifold_id_,
                     int                        levels_,
                     double                     distance_max_)
    : GridRefinement(pcout)
    , _box(box)
    , manifold_id(manifold_id_)
    , levels(levels_)
    , distance_max(distance_max_){};
  virtual ~DistanceRefinement(){};

  virtual void
  refine(dealii::Triangulation<2, 3> &tria, const std::vector<TopoDS_Shape> &cad_surfaces) override;

private:
  Bnd_Box      _box;
  unsigned int manifold_id  = 0;
  int          levels       = 0;
  double       distance_max = 0.0;
};
#endif // DISTANCE_REFINEMENT_H
