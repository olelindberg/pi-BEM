#ifndef GEOMETRY_INPUT_H
#define GEOMETRY_INPUT_H

#include <string>

class GeometryInput
{
public:
  void
  set_filename(std::string val)
  {
    _filename = val;
  }

  void
  set_scale(double val)
  {
    _scale = val;
  }

  void
  set_rotate_x(double val)
  {
    _rotate_x = val;
  }

  void
  set_rotate_y(double val)
  {
    _rotate_y = val;
  }

  void
  set_rotate_z(double val)
  {
    _rotate_z = val;
  }

  void
  set_pre_translate_x(double val)
  {
    _pre_translate_x = val;
  }

  void
  set_pre_translate_y(double val)
  {
    _pre_translate_y = val;
  }

  void
  set_pre_translate_z(double val)
  {
    _pre_translate_z = val;
  }

  void
  set_post_translate_x(double val)
  {
    _post_translate_x = val;
  }

  void
  set_post_translate_y(double val)
  {
    _post_translate_y = val;
  }

  void
  set_post_translate_z(double val)
  {
    _post_translate_z = val;
  }

  std::string
  get_filename() const
  {
    return _filename;
  }

  double
  get_scale() const
  {
    return _scale;
  }

  double
  get_rotate_x() const
  {
    return _rotate_x;
  }

  double
  get_rotate_y() const
  {
    return _rotate_y;
  }

  double
  get_rotate_z() const
  {
    return _rotate_z;
  }

  double
  get_pre_translate_x() const
  {
    return _pre_translate_x;
  }

  double
  get_pre_translate_y() const
  {
    return _pre_translate_y;
  }

  double
  get_pre_translate_z() const
  {
    return _pre_translate_z;
  }

  double
  get_post_translate_x() const
  {
    return _post_translate_x;
  }

  double
  get_post_translate_y() const
  {
    return _post_translate_y;
  }

  double
  get_post_translate_z() const
  {
    return _post_translate_z;
  }

private:
  std::string _filename = "none";

  double _pre_translate_x = 0.0;
  double _pre_translate_y = 0.0;
  double _pre_translate_z = 0.0;

  double _scale = 1.0;

  double _rotate_x = 0.0;
  double _rotate_y = 0.0;
  double _rotate_z = 0.0;

  double _post_translate_x = 0.0;
  double _post_translate_y = 0.0;
  double _post_translate_z = 0.0;
};

#endif