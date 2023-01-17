#ifndef SHAPE_INPUT_H
#define SHAPE_INPUT_H


#include "GeometryInput.h"

class ShapeInput : public GeometryInput
{
public:
  void
  set_mesh_projection(std::string val)
  {
    _mesh_projection = val;
  }

  void
  set_mesh_element_id(int val)
  {
    _mesh_element_id = val;
  }

  std::string
  get_mesh_projection() const
  {
    return _mesh_projection;
  }

  int
  get_mesh_element_id() const
  {
    return _mesh_element_id;
  }

private:
  std::string _mesh_projection = "none";
  int         _mesh_element_id = 0;
};

#endif