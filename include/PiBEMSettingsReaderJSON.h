#ifndef PIBEM_SETTINGS_READER_JSON_H
#define PIBEM_SETTINGS_READER_JSON_H

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <exception>
#include <iostream>
#include <vector>

#include "PiBEMSettings.h"



class PiBEMSettingsReaderJSON
{
public:
  static bool read(const std::string &filename, PiBEMSettings &pibemSettings)
  {
    namespace pt = boost::property_tree;
    pt::ptree root;
    try
    {
      pt::read_json(filename, root);

      pibemSettings.gravity        = root.get<double>("gravity");
      pibemSettings.density        = root.get<double>("density");
      pibemSettings.aspectRatioMax = root.get<double>("adaptive_mesh_refinement.aspectRatioMax");
      pibemSettings.cellSizeMin    = root.get<double>("adaptive_mesh_refinement.cellSizeMin");
      pibemSettings.adaptiveRefinementLevels =
        root.get<double>("adaptive_mesh_refinement.adaptiveRefinementLevels");
      pibemSettings.iterMax = root.get<double>("adaptive_mesh_refinement.iterMax");
      pibemSettings.number_of_elements_max =
        root.get<int>("adaptive_mesh_refinement.number_of_elements_max");
      pibemSettings.top_fraction_max =
        root.get<double>("adaptive_mesh_refinement.top_fraction_max");

      return true;
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      return false;
    }
  }

private:
};

#endif // PIBEM_SETTINGS_READER_JSON_H
