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
  static bool
  read(const std::string &filename, PiBEMSettings &pibemSettings)
  {
    namespace pt = boost::property_tree;
    pt::ptree root;
    try
    {
      pt::read_json(filename, root);

      pibemSettings.gravity                    = root.get<double>("gravity");
      pibemSettings.density                    = root.get<double>("density");
      pibemSettings.aspectRatioMax             = root.get<double>("aspectRatioMax");
      pibemSettings.potentialErrorEstimatorMax = root.get<double>("potentialErrorEstimatorMax");
      pibemSettings.velocityErrorEstimatorMax  = root.get<double>("velocityErrorEstimatorMax");
      pibemSettings.adaptiveRefinementLevels   = root.get<double>("adaptiveRefinementLevels");

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
