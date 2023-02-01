#ifndef BVH_SURFACESELECTOR_H
#define BVH_SURFACESELECTOR_H

#include <BRep_Tool.hxx>
#include <BVH_BoxSet.hxx>
#include <BVH_Traverse.hxx>
#include <GC_MakeLine.hxx>
#include <GeomAPI_IntCS.hxx>
#include <Geom_Line.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Lin.hxx>



typedef BVH_Box<double, 2>                 bvh_box_t;
typedef BVH_BoxSet<double, 2, TopoDS_Face> bvh_boxset_t;
typedef bvh_box_t::BVH_VecNt               bvh_vec_t;



class BVH_SurfaceSelector : public BVH_Traverse<double, 2, bvh_boxset_t, double>
{
public:
  virtual Standard_Boolean
  RejectNode(const bvh_vec_t &theCMin, const bvh_vec_t &theCMax, double &) const Standard_OVERRIDE
  {
    if (theCMin.x() < _point.x() && _point.x() < theCMax.x())
      if (theCMin.y() < _point.y() && _point.y() < theCMax.y())
        return false;
    return true;
  }

  virtual Standard_Boolean Accept(const Standard_Integer theIndex, const double &) Standard_OVERRIDE
  {
    const TopoDS_Face &face = myBVHSet->Element(theIndex);
    const bvh_box_t &  box  = myBVHSet->Box(theIndex);

    double dummy;
    if (!this->RejectNode(box.CornerMin(), box.CornerMax(), dummy))
    {
      _surfaces.push_back(BRep_Tool::Surface(myBVHSet->Element(theIndex)));
      return true;
    }
    return false;
  }

  void set_point(const bvh_vec_t &point)
  {
    _point = point;
  }
  const std::vector<Handle(Geom_Surface)> &get_surfaces()
  {
    return _surfaces;
  }

private:
  bvh_vec_t                         _point;
  std::vector<Handle(Geom_Surface)> _surfaces;
};
#endif // BVH_SURFACESELECTOR_H