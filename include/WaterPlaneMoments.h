#ifndef WATER_PLANE_MOMENTS_H
#define WATER_PLANE_MOMENTS_H

class WaterPlaneMoments
{
public:
  double
  getS0() const
  {
    return _S0;
  }
  double
  getSx() const
  {
    return _Sx;
  }
  double
  getSy() const
  {
    return _Sy;
  }
  double
  getSxx() const
  {
    return _Sxx;
  }
  double
  getSxy() const
  {
    return _Sxy;
  }
  double
  getSyy() const
  {
    return _Syy;
  }

  void
  setS0(double val)
  {
    _S0 = val;
  }
  void
  setSx(double val)
  {
    _Sx = val;
  }
  void
  setSy(double val)
  {
    _Sy = val;
  }
  void
  setSxx(double val)
  {
    _Sxx = val;
  }
  void
  setSxy(double val)
  {
    _Sxy = val;
  }
  void
  setSyy(double val)
  {
    _Syy = val;
  }

private:
  double _S0  = 0.0;
  double _Sx  = 0.0;
  double _Sy  = 0.0;
  double _Sxx = 0.0;
  double _Sxy = 0.0;
  double _Syy = 0.0;
};

#endif // WATER_PLANE_MOMENTS_H
