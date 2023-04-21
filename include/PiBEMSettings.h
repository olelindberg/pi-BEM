#ifndef PIBEM_SETTINGS_H
#define PIBEM_SETTINGS_H

#include <iostream>

class PiBEMSettings
{
public:
  enum class DOMAIN_TYPE
  {
    SINGLE_MESH,
    MULTI_MESH
  };

  double      gravity                  = 9.80665;
  double      density                  = 1000;
  double      aspectRatioMax           = 2.5;
  double      cellSizeMin              = 0.0;
  int         adaptiveRefinementLevels = 0;
  int         iterMax                  = 0;
  int         number_of_elements_max   = 2000;
  double      top_fraction_max         = 0.05;
  DOMAIN_TYPE domain_type              = DOMAIN_TYPE::SINGLE_MESH;
  bool        route_enabled            = true;

  void print()
  {
    std::cout << "\npi-BEM settings:" << std::endl;
    std::cout << "gravity                   = " << gravity << std::endl;
    std::cout << "density                   = " << density << std::endl;
    std::cout << "aspectRatioMax            = " << aspectRatioMax << std::endl;
    std::cout << "cellSizeMin               = " << cellSizeMin << std::endl;
    std::cout << "adaptiveRefinementLevels  = " << adaptiveRefinementLevels << std::endl;
    std::cout << "iterMax                   = " << iterMax << std::endl;
    std::cout << "number_of_elements_max    = " << number_of_elements_max << std::endl;
    std::cout << "top_fraction_max          = " << top_fraction_max << std::endl;
    std::cout << "domain_type               = " << (int)domain_type << std::endl;
    std::cout << "route_enabled             = " << (int)route_enabled << std::endl;
  }
};

#endif // PIBEM_SETTINGS_H
