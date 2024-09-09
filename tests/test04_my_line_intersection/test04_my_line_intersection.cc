#include "my_utilities.h"
#include <deal.II/opencascade/utilities.h>
#include <iostream>

#include <ShapeAnalysis_Surface.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
//#include <Bnd_OBB.hxx>
//#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>

int
main ()
{
  //dealii::Point<3>     origin (-0.642229, -0.642229, -0.288675);
  //dealii::Tensor<1, 3> direct;
  //direct[0] = 0.672922;
  //direct[1] = 0.674851;
  //direct[2] = 0.302905;
//
  //std::string      filename  = "/home/ole/dev/projects/pi-BEM/examples/example04_hemisphere/Color_1.iges";
  //TopoDS_Shape     shape     = dealii::OpenCASCADE::read_IGES (filename, 1e-3);
  //double           tolerance = 1.0e-13;
  //dealii::Point<3> point     = my_line_intersection (shape, origin, direct, tolerance);
//
  //std::cout << "Intersection point " << point << std::endl;
  //{
//
  //  gp_Pnt              P (origin[0], origin[1], origin[2]);
  //  std::vector<gp_Pnt> surfacePoints1;
  //  std::vector<gp_Pnt> surfacePoints2;
  //  std::vector<double> dist1;
  //  std::vector<double> dist2;
  //  TopExp_Explorer     Ex;
  //  for (Ex.Init (shape, TopAbs_FACE); Ex.More (); Ex.Next ())
  //  {
  //    Handle (Geom_Surface) surface = BRep_Tool::Surface (TopoDS::Face (Ex.Current ()));
//
  //    ShapeAnalysis_Surface projector (surface);
  //    gp_Pnt2d              surfaceParam = projector.ValueOfUV (P, 1.0e-10);
//
  //    gp_Pnt surfacePoint;
  //    surface->D0 (surfaceParam.X (), surfaceParam.Y (), surfacePoint);
  //    dist1.push_back( P.Distance (surfacePoint));
  //    surfacePoints1.push_back (surfacePoint);
//
  //    GeomAPI_ProjectPointOnSurf projector2 (P, surface);
  //    for (int i = 0; i < projector2.NbPoints (); ++i)
  //    {
  //    dist2.push_back( P.Distance (projector2.Point(i + 1)));
  //      surfacePoints2.push_back (projector2.Point(i + 1));
  //    }
  //  }
//
  //  for (auto& sp : surfacePoints1)
  //    std::cout << "surfacePoints1 " << sp.X() << ", " << sp.Y() << ", " << sp.Z()  << std::endl;
  //  for (auto& d : dist1)
  //    std::cout << "dist1 " << d  << std::endl;
//
  //  for (auto& sp : surfacePoints2)
  //    std::cout << "surfacePoints2 " << sp.X() << ", " << sp.Y() << ", " << sp.Z()  << std::endl;
  //  for (auto& d : dist2)
  //    std::cout << "dist2 " << d  << std::endl;
  //}
  return 0;
}
