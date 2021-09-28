#ifndef GRID_REFINEMENT_CREATOR_H
#define GRID_REFINEMENT_CREATOR_H

#include "../include/IGridRefinement.h"

#include <deal.II/base/conditional_ostream.h>

#include <TopoDS_Shape.hxx>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>


#include <exception>
#include <iostream>
#include <vector>

#include "Body.h"

class GridRefinementCreator
{
public:
  static std::vector<std::shared_ptr<IGridRefinement>>
  create(const std::string &              filename,
         dealii::ConditionalOStream       pcout,
         double                           tolerance,
         const std::vector<TopoDS_Shape> &cad_surfaces);

private:
};

#endif // GRID_REFINEMENT_CREATOR_H
