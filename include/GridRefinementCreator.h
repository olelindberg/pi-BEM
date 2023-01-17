#ifndef GRID_REFINEMENT_CREATOR_H
#define GRID_REFINEMENT_CREATOR_H

#include "../include/IGridRefinement.h"

#include <deal.II/base/conditional_ostream.h>

#include <TopoDS_Shape.hxx>



#include <exception>
#include <iostream>
#include <vector>

#include "Body.h"

class GridRefinementCreator
{
public:
  static std::vector<std::shared_ptr<IGridRefinement>>
  create(const std::string &filename, dealii::ConditionalOStream pcout, double tolerance);

private:
};

#endif // GRID_REFINEMENT_CREATOR_H
