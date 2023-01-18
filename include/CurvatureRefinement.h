#ifndef CURVATURE_REFINEMENT_H
#define CURVATURE_REFINEMENT_H

#include "GridRefinement.h"

#include <iostream>
#include <vector>

class CurvatureRefinement : public GridRefinement
{
public:
  CurvatureRefinement(dealii::ConditionalOStream pcout,
                      bool                       use_cad_surface_and_curves_,
                      bool                       surface_curvature_refinement_,
                      double                     cells_per_circle_,
                      double                     cad_to_projectors_tolerance_ratio_,
                      double                     max_tol_,
                      unsigned int               max_curvature_ref_cycles_)
    : GridRefinement(pcout)
    , use_cad_surface_and_curves(use_cad_surface_and_curves_)
    , surface_curvature_refinement(surface_curvature_refinement_)
    , cells_per_circle(cells_per_circle_)
    , cad_to_projectors_tolerance_ratio(cad_to_projectors_tolerance_ratio_)
    , max_tol(max_tol_)
    , max_curvature_ref_cycles(max_curvature_ref_cycles_)



  {}

  virtual ~CurvatureRefinement()
  {}

  virtual void refine(dealii::Triangulation<2, 3> &tria) override;

private:
  bool         use_cad_surface_and_curves        = true;
  bool         surface_curvature_refinement      = true;
  double       cells_per_circle                  = 0.0;
  double       cad_to_projectors_tolerance_ratio = 0.0;
  double       max_tol                           = 0.0;
  unsigned int max_curvature_ref_cycles          = 0;
};

#endif // CURVATURE_REFINEMENT_H