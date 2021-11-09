#include "../include/GridRefinementCreator.h"

#include "../include/AspectRatioRefinement.h"
#include "../include/CurvatureRefinement.h"
#include "../include/DistanceRatioRefinement.h"
#include "../include/HeightRatioRefinement.h"
#include "../include/ManifoldRefinement.h"

#include <BRepBndLib.hxx>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

std::vector<std::shared_ptr<IGridRefinement>>
GridRefinementCreator::create(const std::string &              filename,
                              dealii::ConditionalOStream       pcout,
                              double                           tolerance,
                              const std::vector<TopoDS_Shape> &cad_surfaces)
{
  std::vector<std::shared_ptr<IGridRefinement>> gridrefinement;

  namespace pt = boost::property_tree;
  pt::ptree root;
  try
  {
    pt::read_json(filename, root);
    {
      auto child = root.get_child("aspectRatioRefinement");
      if (!child.empty())
      {
        pcout << "Creating aspect ratio refinement ...\n";
        double           aspect_ratio_max = child.get<double>("aspect_ratio_max");
        std::vector<int> manifold_ids;
        for (pt::ptree::value_type &id : child.get_child("manifold_ids"))
          manifold_ids.push_back(id.second.get_value<int>());
        double refinement_levels = child.get<double>("refinement_levels");
        gridrefinement.push_back(std::make_shared<AspectRatioRefinement>(
          pcout, manifold_ids, refinement_levels, aspect_ratio_max));
      }
    }

    {
      auto child = root.get_child("heightRatioRefinement");
      if (!child.empty())
      {
        pcout << "Creating height ratio refinement ...\n";
        double height_ratio_max  = child.get<double>("height_ratio_max");
        double refinement_levels = child.get<double>("refinement_levels");
        gridrefinement.push_back(std::make_shared<HeightRatioRefinement>(
          pcout, refinement_levels, height_ratio_max, tolerance));
      }
    }

    {
      auto child = root.get_child("distanceRefinement");
      if (!child.empty())
      {
        pcout << "Creating distance refinement ...\n";
        Bnd_Box box;
        for (const auto &shape : cad_surfaces)
          BRepBndLib::Add(shape, box);
        int    manifold_id       = child.get<int>("manifold_id");
        int    refinement_levels = child.get<int>("refinement_levels");
        double distance_max      = child.get<double>("distance_max");
        gridrefinement.push_back(std::make_shared<DistanceRatioRefinement>(
          pcout, box, manifold_id, refinement_levels, distance_max));
      }
    }
  }
  catch (std::exception &e)
  {
    pcout << e.what() << std::endl;
  }
  return gridrefinement;
}
