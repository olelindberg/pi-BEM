#ifndef WIRE_UTIL_H
#define WIRE_UTIL_H

#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>

#include "SurfaceMoments.h"

class WireUtil
{
public:
  static bool
  surfaceMoments(const TopoDS_Wire &wire, SurfaceMoments &sm)
  {
    BRepBuilderAPI_MakeFace facebuilder(wire);
    if (facebuilder.IsDone())
    {
      GProp_GProps gprops(gp_Pnt(sm.getx0(), sm.gety0(), 0.0));
      BRepGProp::SurfaceProperties(facebuilder.Face(), gprops);
      sm.setS0(gprops.Mass());

      double Sx, Sy, Sz;
      gprops.StaticMoments(Sx, Sy, Sz);
      sm.setSx(Sx);
      sm.setSy(Sy);

      gp_Mat mat = gprops.MatrixOfInertia();
      sm.setSxx(mat.Value(1, 1));
      sm.setSxy(mat.Value(1, 2));
      sm.setSyy(mat.Value(2, 2));
      return true;
    }
    return false;
  }
};
#endif // WIRE_UTIL_H