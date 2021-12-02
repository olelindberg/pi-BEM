#ifndef PIBEM_SETTINGS_H
#define PIBEM_SETTINGS_H

#include <iostream>

class PiBEMSettings
{
public:
  double gravity                    = 9.80665;
  double density                    = 1000;
  double potentialErrorEstimatorMax = 1.25;
  double velocityErrorEstimatorMax  = 1.25;
  double aspectRatioMax             = 2.5;
  double cellSizeMin                = 0.0;
  int    adaptiveRefinementLevels   = 0;
  int    iterMax                    = 0;

  void
  print()
  {
    std::cout << "pi-BEM settings:" << std::endl;
    std::cout << "gravity                     = " << gravity << std::endl;
    std::cout << "density                     = " << density << std::endl;
    std::cout << "potentialErrorEstimatorMax  = " << potentialErrorEstimatorMax << std::endl;
    std::cout << "velocityErrorEstimatorMax   = " << velocityErrorEstimatorMax << std::endl;
    std::cout << "aspectRatioMax              = " << aspectRatioMax << std::endl;
    std::cout << "cellSizeMin                 = " << cellSizeMin << std::endl;
    std::cout << "adaptiveRefinementLevels    = " << adaptiveRefinementLevels << std::endl;
    std::cout << "iterMax                     = " << iterMax << std::endl;
  }
};

#endif // PIBEM_SETTINGS_H
