#ifndef SURFACE_MOMENTS_H
#define SURFACE_MOMENTS_H

class SurfaceMoments
{
public:
  SurfaceMoments(double x0, double y0)
    : _x0(x0)
    , _y0(y0)
  {}

  double getx0() const
  {
    return _x0;
  }
  double gety0() const
  {
    return _y0;
  }
  double getS0() const
  {
    return _S0;
  }

  double getSx() const
  {
    return _Sx;
  }
  double getSy() const
  {
    return _Sy;
  }
  double getSxx() const
  {
    return _Sxx;
  }
  double getSxy() const
  {
    return _Sxy;
  }
  double getSyy() const
  {
    return _Syy;
  }

  void setS0(double val)
  {
    _S0 = val;
  }
  void setSx(double val)
  {
    _Sx = val;
  }
  void setSy(double val)
  {
    _Sy = val;
  }
  void setSxx(double val)
  {
    _Sxx = val;
  }
  void setSxy(double val)
  {
    _Sxy = val;
  }
  void setSyy(double val)
  {
    _Syy = val;
  }

  void print()
  {
    std::cout << "x0  : " << _x0 << std::endl;
    std::cout << "y0  : " << _y0 << std::endl;
    std::cout << "S0  : " << _S0 << std::endl;
    std::cout << "Sx  : " << _Sx << std::endl;
    std::cout << "Sy  : " << _Sy << std::endl;
    std::cout << "Sxx : " << _Sxx << std::endl;
    std::cout << "Sxy : " << _Sxy << std::endl;
    std::cout << "Syy : " << _Syy << std::endl;
  }

private:
  double _x0  = 0.0;
  double _y0  = 0.0;
  double _S0  = 0.0;
  double _Sx  = 0.0;
  double _Sy  = 0.0;
  double _Sxx = 0.0;
  double _Sxy = 0.0;
  double _Syy = 0.0;
};

#endif // SURFACE_MOMENTS_H
