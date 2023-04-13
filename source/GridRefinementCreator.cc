#include "../include/GridRefinementCreator.h"

#include "../include/ErrorMessage.h"

#include "../include/AreaRefinement.h"
#include "../include/AspectRatioRefinement.h"
#include "../include/BoxRefinement.h"
#include "../include/CurvatureRefinement.h"
#include "../include/DistanceRatioRefinement.h"
#include "../include/HeightRatioRefinement.h"
#include "../include/ManifoldRefinement.h"

#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <gp_Pnt.hxx>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

std::vector<std::shared_ptr<IGridRefinement>>
GridRefinementCreator::create(const std::string &filename, dealii::ConditionalOStream pcout)
{
  std::vector<std::shared_ptr<IGridRefinement>> gridrefinement;

  namespace pt = boost::property_tree;
  pt::ptree root;
  try
  {
    pt::read_json(filename, root);
  }
  catch (std::exception &e)
  {
    std::cout << ErrorMessage::message(__FILE__, __LINE__, e.what());
    return gridrefinement;
  }

  try
  {
    for (auto &node : root)
    {
      if (node.first == "aspectRatioRefinement")
      {
        auto child = node.second;
        if (!child.empty())
        {
          double           aspect_ratio_max = child.get<double>("aspect_ratio_max");
          std::vector<int> manifold_ids;
          for (pt::ptree::value_type &id : child.get_child("manifold_ids"))
            manifold_ids.push_back(id.second.get_value<int>());
          double refinement_levels = child.get<double>("refinement_levels");
          bool   verbose           = child.get<bool>("verbose", false);
          gridrefinement.push_back(std::make_shared<AspectRatioRefinement>(
            pcout, manifold_ids, refinement_levels, aspect_ratio_max, verbose));
        }
      }

      if (node.first == "heightRatioRefinement")
      {
        auto child = node.second;
        if (!child.empty())
        {
          double height_ratio_max  = child.get<double>("height_ratio_max");
          double refinement_levels = child.get<double>("refinement_levels");
          gridrefinement.push_back(
            std::make_shared<HeightRatioRefinement>(pcout, refinement_levels, height_ratio_max));
        }
      }

      if (node.first == "areaRefinement")
      {
        auto child = node.second;
        if (!child.empty())
        {
          std::vector<int> manifold_ids;
          for (pt::ptree::value_type &id : child.get_child("manifold_ids"))
            manifold_ids.push_back(id.second.get_value<int>());

          int    refinement_levels = child.get<int>("refinement_levels");
          double area_max          = child.get<double>("area_max");
          double aspect_ratio_max  = child.get<double>("aspect_ratio_max");
          bool   verbose           = child.get<bool>("verbose", false);

          gridrefinement.push_back(std::make_shared<AreaRefinement>(
            pcout, manifold_ids, aspect_ratio_max, refinement_levels, area_max, verbose));
        }
      }


      if (node.first == "distanceRatioRefinement")
      {
        auto child = node.second;
        if (!child.empty())
        {
          int    manifold_id        = child.get<int>("manifold_id");
          int    refinement_levels  = child.get<int>("refinement_levels");
          double distance_ratio_max = child.get<double>("distance_ratio_max");
          double aspect_ratio_max   = child.get<double>("aspect_ratio_max");
          double xmin               = child.get<double>("xmin");
          double xmax               = child.get<double>("xmax");
          double ymin               = child.get<double>("ymin");
          double ymax               = child.get<double>("ymax");
          double zmin               = child.get<double>("zmin");
          double zmax               = child.get<double>("zmax");

          gp_Pnt pmin(xmin, ymin, zmin);
          gp_Pnt pmax(xmax, ymax, zmax);

          Bnd_Box box;
          box.Add(pmin);
          box.Add(pmax);


          gridrefinement.push_back(std::make_shared<DistanceRatioRefinement>(
            pcout, box, aspect_ratio_max, manifold_id, refinement_levels, distance_ratio_max));
        }
      }

      if (node.first == "boxRefinement")
      {
        auto child = node.second;
        if (!child.empty())
        {
          int    manifold_id       = child.get<int>("manifold_id");
          int    refinement_levels = child.get<int>("refinement_levels");
          double aspect_ratio_max  = child.get<double>("aspect_ratio_max");
          double area_min          = child.get<double>("area_min");

          double xmin = child.get<double>("xmin");
          double xmax = child.get<double>("xmax");
          double ymin = child.get<double>("ymin");
          double ymax = child.get<double>("ymax");
          double zmin = child.get<double>("zmin");
          double zmax = child.get<double>("zmax");

          gp_Pnt pmin(xmin, ymin, zmin);
          gp_Pnt pmax(xmax, ymax, zmax);

          Bnd_Box box;
          box.Add(pmin);
          box.Add(pmax);

          gridrefinement.push_back(std::make_shared<BoxRefinement>(
            pcout, box, aspect_ratio_max, area_min, manifold_id, refinement_levels));
        }
      }
    }
  }
  catch (const std::exception &e)
  {
    std::cout << ErrorMessage::message(__FILE__, __LINE__, e.what());
  }

  return gridrefinement;
}
