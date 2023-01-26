#ifndef KCS_LIMFJORD_SETUP_H
#define KCS_LIMFJORD_SETUP_H

#include "GeometrySetup.h"
#include <filesystem>

class KCS_LimfjordSetup : public GeometrySetup
{
public:
  KCS_LimfjordSetup()
  {
    //---------------------------------------------------------------------------
    // Mesh input:
    //---------------------------------------------------------------------------
    {
      GeometryInput input;
      input.set_filename("mesh1.inp");
      _mesh_inputs.push_back(input);
    }

    {
      GeometryInput input;
      input.set_filename("mesh2.inp");
      _mesh_inputs.push_back(input);
    }

    //---------------------------------------------------------------------------
    // CAD geometry input:
    //---------------------------------------------------------------------------

    {
      ShapeInput input;
      input.set_filename("Color_1.iges");
      input.set_mesh_projection("normal");
      input.set_mesh_element_id(1);
      _shape_inputs.push_back(input);
    }

    {
      ShapeInput input;
      input.set_filename("Color_2.iges");
      input.set_mesh_projection("normal");
      input.set_mesh_element_id(2);
      _shape_inputs.push_back(input);
    }

    {
      ShapeInput input;
      input.set_filename("limfjord_seabed.iges");
      input.set_mesh_projection("directional");
      input.set_mesh_element_id(3);
      //      input.set_pre_translate_x(-581628.0);
      //      input.set_pre_translate_y(-6314958.0);
      input.set_pre_translate_z(10.8);
      //      input.set_rotate_z(-2.5);

      _shape_inputs.push_back(input);
    }

    {
      ShapeInput input;
      input.set_filename("Curve_1.iges");
      input.set_mesh_projection("arch_length");
      input.set_mesh_element_id(11);
      _shape_inputs.push_back(input);
    }

    {
      ShapeInput input;
      input.set_filename("Curve_2.iges");
      input.set_mesh_projection("arch_length");
      input.set_mesh_element_id(12);
      _shape_inputs.push_back(input);
    }

    {
      ShapeInput input;
      input.set_filename("Curve_3.iges");
      input.set_mesh_projection("arch_length");
      input.set_mesh_element_id(13);
      _shape_inputs.push_back(input);
    }

    {
      ShapeInput input;
      input.set_filename("Curve_4.iges");
      input.set_mesh_projection("arch_length");
      input.set_mesh_element_id(14);
      _shape_inputs.push_back(input);
    }
  }


  virtual std::vector<GeometryInput> mesh_inputs() override
  {
    return _mesh_inputs;
  }

  virtual std::vector<ShapeInput> shape_inputs() override
  {
    return _shape_inputs;
  }

private:
  std::vector<ShapeInput>    _shape_inputs;
  std::vector<GeometryInput> _mesh_inputs;
};

#endif // KCS_LIMFJORD_SETUP_H
