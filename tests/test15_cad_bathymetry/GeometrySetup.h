#ifndef GEOMETRY_SETUP_H
#define GEOMETRY_SETUP_H

#include "ShapeInput.h"
#include <vector>

class GeometrySetup
{
public:
  virtual std::vector<GeometryInput>
  mesh_inputs() = 0;
  virtual std::vector<ShapeInput>
  shape_inputs() = 0;
};

#endif // GEOMETRY_SETUP_H
